/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2026  Laurent Forthomme
 *                2017-2019  Wolfgang Schaefer
 *                2019       Marta Luszczak
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "PPtoPair.h"

using namespace cepgen;
using namespace cepgen::proc;
using namespace std::string_literals;

PPtoPair::PPtoPair(const ParametersList& params, const std::vector<short>& hel_x1, const std::vector<short>& hel_x2)
    : cepgen::proc::FactorisedProcess(params, buildPdgIdsPair(params.get<ParticleProperties>("pair").pdgid)),
      osp_(steer<ParametersList>("offShellParameters")),
      helicities_x1_(hel_x1),
      helicities_x2_(hel_x2),
      method_(steerAs<int, Mode>("method")) {}

ParametersDescription PPtoPair::description() {
  auto desc = FactorisedProcess::description();
  desc.setDescription("γγ → XXbar");
  desc.addAs<int, pdgid_t>("pair", PDG::muon).setDescription("type of central particles emitted");
  desc.addAs<int, Mode>("method", Mode::offShell)
      .setDescription("Matrix element computation method")
      .allow(0, "on-shell")
      .allow(1, "off-shell")
      .allow(2, "off-shell (helicity amplitudes)");
  desc.add("offShellParameters", OffShellParameters::description());
  desc.add("model", amplitudes::HelicityAmplitudes::description())
      .setDescription("Helicity amplitudes model (only used for the method 2)");
  return desc;
}

void PPtoPair::prepareFactorisedPhaseSpace() {
  if ((method_ == Mode::offShell || method_ == Mode::helicity) && !partonsKtFactorised())
    throw CG_FATAL("PPtoFF:prepare") << "Off-shell matrix element only defined for factorised process with partons kt.";
  const auto cs_prop = PDG::get()(central().at(0));
  // define central particle properties and couplings with partons
  if (method_ == Mode::helicity)
    amplitudes_ = amplitudes::HelicityAmplitudesFactory::get().build(steer<ParametersList>("model"));
  mx2_ = cs_prop.mass * cs_prop.mass;
  qx2_ = std::pow(cs_prop.integerCharge() * (1. / 3), 2);
  const auto generate_coupling = [this, &cs_prop](const pdgid_t& parton_id) -> std::function<double(double)> {
    switch (parton_id) {
      case PDG::gluon: {
        if (cs_prop.colours == 0)
          throw CG_FATAL("PPtoFF:prepare") << "Invalid fermion type for gluon coupling. Should be a quark.";
        return [this](double q) { return kFourPi * 0.5 * alphaS(q); };
      }
      case PDG::photon:
        return [this](double q) { return kFourPi * qx2_ * alphaEM(q); };
      default:
        throw CG_FATAL("PPtoFF:prepare") << "Unsupported parton id: '" << parton_id << "'.";
    }
  };
  const auto parton1_id = event().oneWithRole(Particle::Role::Parton1).pdgId(),
             parton2_id = event().oneWithRole(Particle::Role::Parton2).pdgId();
  g_part1_ = generate_coupling(parton1_id);
  g_part2_ = generate_coupling(parton2_id);
  // sum over final-state colours for colour-neutral (photon-photon) initial states
  // (colour factor already accounted for in coupling for gluon-induced processes)
  colour_factor_ =
      (parton1_id == PDG::photon && parton2_id == PDG::photon && cs_prop.colours > 0) ? cs_prop.colours : 1.;
}

double PPtoPair::offShellME() const {
  if (q1().pt2() == 0. || q2().pt2() == 0)  // only works for kt-factorised case
    return 0.;
  const auto mt1 = pc(0).massT(), mt2 = pc(1).massT();  // transverse masses
  const auto compute_zs = [this, &mt1, &mt2](short pol, double x) -> std::pair<double, double> {
    const auto norm_pol = pol / std::abs(pol);
    const auto fact = inverseSqrtS() / x;
    return std::make_pair(fact * mt1 * std::exp(norm_pol * pc(0).rapidity()),
                          fact * mt2 * std::exp(norm_pol * pc(1).rapidity()));
  };
  //--- t-channel
  const auto q2_1 = utils::kt::q2(x1(), q1().pt2(), mA2(), mX2());
  const auto [zp_1, zm_1] = compute_zs(+1, x1());
  const auto amat2_1 = phiMatrixElement(zp_1, zm_1, q2_1, q1(), q2().transverse());

  //--- u-channel
  const auto q2_2 = utils::kt::q2(x2(), q2().pt2(), mB2(), mY2());
  const auto [zp_2, zm_2] = compute_zs(-1, x2());
  const auto amat2_2 = phiMatrixElement(zp_2, zm_2, q2_2, q2(), q1().transverse());

  //--- symmetrisation
  const auto amat2 = 0.5 * (osp_.mat1 * amat2_1 + osp_.mat2 * amat2_2);
  if (!utils::positive(amat2))
    return 0.;

  const auto t_limits = Limits{0., std::pow(std::max(mt1, mt2), 2)};
  if (const auto prefactor =
          colour_factor_ * g_part1_(std::sqrt(t_limits.trim(q2_1))) * g_part2_(std::sqrt(t_limits.trim(q2_2)));
      utils::positive(prefactor))
    return prefactor * amat2;
  return 0.;
}

double PPtoPair::helicityME() const {
  // off-shell matrix element built from the on-shell γγ → f f̄ helicity amplitudes, with the incoming photons
  // polarised along their respective transverse momenta (see e.g. \cite Luszczak:2018ntp for the γγ → W⁺W⁻ case)
  if (q1().pt2() == 0. || q2().pt2() == 0.)
    return 0.;
  const auto s_hat = shat();
  if (s_hat <= 4. * mx2_)
    return 0.;
  const amplitudes::Kinematics kin(mx2_, s_hat, that(), uhat());
  // the helicity amplitudes are defined in the scattering plane: polarisations are taken relative to the pair azimuth
  const auto phi_pair = std::atan2(pc(0).py() - pc(1).py(), pc(0).px() - pc(1).px());
  const auto e1 = std::polar(1., std::atan2(q1().py(), q1().px()) - phi_pair),
             e2 = std::polar(1., std::atan2(q2().py(), q2().px()) - phi_pair);
  double amat2{0.};
  for (const short lam3 : helicities_x1_)
    for (const short lam4 : helicities_x2_)
      amat2 += std::norm(amplitudes_->amplitudeLinearPolarisations(kin, lam3, lam4, e1, e2));
  amat2 /= s_hat * s_hat;
  if (!utils::positive(amat2))
    return 0.;

  const auto mt1 = pc(0).massT(), mt2 = pc(1).massT();  // transverse masses
  const auto q2_1 = utils::kt::q2(x1(), q1().pt2(), mA2(), mX2()), q2_2 = utils::kt::q2(x2(), q2().pt2(), mB2(), mY2());
  const auto t_limits = Limits{0., std::pow(std::max(mt1, mt2), 2)};
  if (const auto prefactor = couplingFactor(std::sqrt(t_limits.trim(q2_1)), std::sqrt(t_limits.trim(q2_2)));
      utils::positive(prefactor))
    return prefactor * amat2;
  return 0.;
}

double PPtoPair::couplingFactor(double q1, double q2) const {
  return colour_factor_ * kFourPi * kFourPi * qx2_ * qx2_ * alphaEM(q1) * alphaEM(q2);
}

/// Parameters for the off-shell matrix element
PPtoPair::OffShellParameters::OffShellParameters(const ParametersList& params) : SteeredObject(params) {
  (*this)
      .add("mat1"s, mat1)
      .add("mat2"s, mat2)
      .add("termLL"s, term_ll)
      .add("termLT"s, term_lt)
      .add("termTT"s, term_tt1)
      .add("termtt"s, term_tt2);
}

ParametersDescription PPtoPair::OffShellParameters::description() {
  auto desc = ParametersDescription();
  desc.add("mat1"s, 1).setDescription("symmetry factor for the first incoming photon");
  desc.add("mat2"s, 1).setDescription("symmetry factor for the second incoming photon");
  desc.add("termLL"s, 1).setDescription("fully longitudinal relative weight");
  desc.add("termLT"s, 1).setDescription("cross-polarisation relative weight");
  desc.add("termTT"s, 1).setDescription("fully transverse relative weight");
  desc.add("termtt"s, 1).setDescription("transverse-transverse azimuthal correlation relative weight");
  return desc;
}
