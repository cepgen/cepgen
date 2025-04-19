/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2025  Laurent Forthomme
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
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;

auto make_pdgids_pair = [](pdgid_t pair) { return spdgids_t{(spdgid_t)pair, -(spdgid_t)pair}; };

/// Compute the 2-to-4 matrix element for a CE \f$\gamma\gamma\rightarrow f\bar f\f$ process
class PPtoFF final : public cepgen::proc::FactorisedProcess {
public:
  explicit PPtoFF(const ParametersList& params)
      : cepgen::proc::FactorisedProcess(params, make_pdgids_pair(params.get<ParticleProperties>("pair").pdgid)),
        method_(steerAs<int, Mode>("method")),
        osp_(steer<ParametersList>("offShellParameters")) {
    if (method_ == Mode::offShell && !phase_space_generator_->ktFactorised())
      throw CG_FATAL("PPtoFF:prepare")
          << "Off-shell matrix element only defined for factorised process with partons kt.";
  }

  proc::ProcessPtr clone() const override { return std::make_unique<PPtoFF>(*this); }

  static ParametersDescription description() {
    auto desc = FactorisedProcess::description();
    desc.setDescription("γγ → f⁺f¯");
    desc.addAs<int, pdgid_t>("pair", PDG::muon).setDescription("type of central particles emitted");
    desc.addAs<int, Mode>("method", Mode::offShell)
        .setDescription("Matrix element computation method")
        .allow(0, "on-shell")
        .allow(1, "off-shell");
    desc.add("offShellParameters", OffShellParameters::description());
    return desc;
  }

private:
  void prepareFactorisedPhaseSpace() override {
    const auto cs_prop = PDG::get()(phase_space_generator_->central().at(0));
    // define central particle properties and couplings with partons
    if (!cs_prop.fermion || cs_prop.charges.empty())
      throw CG_FATAL("PPtoFF:prepare") << "Invalid fermion pair selected: " << cs_prop << ".";
    mf2_ = cs_prop.mass * cs_prop.mass;
    qf2_ = std::pow(cs_prop.integerCharge() * (1. / 3), 2);
    const auto generate_coupling = [this, &cs_prop](const pdgid_t& parton_id) -> std::function<double(double)> {
      switch (parton_id) {
        case PDG::gluon: {
          if (cs_prop.colours == 0)
            throw CG_FATAL("PPtoFF:prepare") << "Invalid fermion type for gluon coupling. Should be a quark.";
          return [this](double q) { return kFourPi * 0.5 * alphaS(q); };
        }
        case PDG::photon:
          return [this](double q) { return kFourPi * qf2_ * alphaEM(q); };
        default:
          throw CG_FATAL("PPtoFF:prepare") << "Unsupported parton id: '" << parton_id << "'.";
      }
    };
    g_part1_ = generate_coupling(event().oneWithRole(Particle::Role::Parton1).pdgId());
    g_part2_ = generate_coupling(event().oneWithRole(Particle::Role::Parton2).pdgId());

    CG_DEBUG("PPtoFF:prepare") << "Incoming beams: mA = " << mA() << " GeV/mB = " << mB() << " GeV.\n\t"
                               << "Produced particles: " << phase_space_generator_->central() << ".\n\t"
                               << "ME computation method: " << (int)method_ << ".";

    // constrain central particles cuts
    if (!kinematics().cuts().central.pt_diff.valid())
      kinematics().cuts().central.pt_diff = {0., 50.};
  }
  double computeFactorisedMatrixElement() override {
    switch (method_) {
      case Mode::onShell:
        return onShellME();
      case Mode::offShell:
        return offShellME();
      default:
        throw CG_FATAL("PPtoFF") << "Invalid ME calculation method (" << (int)method_ << ")!";
    }
  }
  double onShellME() const;
  double offShellME() const;

  const enum class Mode { onShell = 0, offShell = 1 } method_;

  /// Parameters for the off-shell matrix element
  struct OffShellParameters : SteeredObject<OffShellParameters> {
    explicit OffShellParameters(const ParametersList& params) : SteeredObject(params) {
      (*this)
          .add("mat1", mat1)
          .add("mat2", mat2)
          .add("termLL", term_ll)
          .add("termLT", term_lt)
          .add("termTT", term_tt1)
          .add("termtt", term_tt2);
    }

    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.add("mat1", 1).setDescription("symmetry factor for the first incoming photon");
      desc.add("mat2", 1).setDescription("symmetry factor for the second incoming photon");
      desc.add("termLL", 1).setDescription("fully longitudinal relative weight");
      desc.add("termLT", 1).setDescription("cross-polarisation relative weight");
      desc.add("termTT", 1).setDescription("fully transverse relative weight");
      desc.add("termtt", 1).setDescription("fully transverse relative weight");
      return desc;
    }

    int mat1{0}, mat2{0};
    int term_ll{0}, term_lt{0}, term_tt1{0}, term_tt2{0};
  } osp_;

  static constexpr double kFourPi = 4. * M_PI;
  double mf2_{0.}, qf2_{0.};
  std::function<double(double)> g_part1_{nullptr}, g_part2_{nullptr};
};

double PPtoFF::onShellME() const {
  const double s_hat = shat(), t_hat = that(), u_hat = uhat();  // buffer Mandelstam variables
  CG_DEBUG_LOOP("PPtoFF:onShell") << "that: " << t_hat << ", uhat: " << u_hat << ".";

  if (t_hat == mf2_ || u_hat == mf2_)
    return 0.;
  const auto q = std::sqrt(t_hat);
  const auto prefac = g_part1_(q) * g_part2_(q);
  if (!utils::positive(prefac))
    return 0.;

  const auto mf4 = mf2_ * mf2_, mf8 = mf4 * mf4;
  const auto out = 6. * mf8 + (-3. * mf4 * t_hat * t_hat) + (-14. * mf4 * t_hat * u_hat) + (-3. * mf4 * u_hat * u_hat) +
                   (1. * mf2_ * t_hat * t_hat * t_hat) + (7. * mf2_ * t_hat * t_hat * u_hat) +
                   (7. * mf2_ * t_hat * u_hat * u_hat) + (1. * mf2_ * u_hat * u_hat * u_hat) +
                   (-1. * t_hat * t_hat * t_hat * u_hat) + (-1. * t_hat * u_hat * u_hat * u_hat);
  return -2. * prefac * out * std::pow((mf2_ - t_hat) * (mf2_ - u_hat) * s_hat, -2);
}

double PPtoFF::offShellME() const {
  if (q1().pt2() == 0. || q2().pt2() == 0)  // only works for kt-factorised case
    return 0.;
  const auto mt1 = pc(0).massT(), mt2 = pc(1).massT();  // transverse masses
  const auto compute_zs = [this, &mt1, &mt2](short pol, double x) -> std::pair<double, double> {
    const auto norm_pol = pol / std::abs(pol);
    const auto fact = inverseSqrtS() / x;
    return std::make_pair(fact * mt1 * std::exp(norm_pol * pc(0).rapidity()),
                          fact * mt2 * std::exp(norm_pol * pc(1).rapidity()));
  };
  const auto compute_mat_element =
      [this](double zp, double zm, double q2, const Momentum& vec_pho, const Momentum& vec_qt) -> double {
    const auto vec_kt = Momentum(zm * pc(0) - zp * pc(1)).transverse();
    const auto phi_p = vec_kt + zp * vec_qt, phi_m = vec_kt - zm * vec_qt;
    const auto zpm = zp * zm, eps2 = mf2_ + zpm * q2;

    const auto kp = 1. / (phi_p.pt2() + eps2), km = 1. / (phi_m.pt2() + eps2);
    const auto phi = Momentum(kp * phi_p - km * phi_m).setEnergy(kp - km);
    const auto dot = phi.threeProduct(vec_pho), cross = phi.crossProduct(vec_pho, Momentum::Z);

    const auto phi_0 = phi.energy(), phi2_0 = phi_0 * phi_0, phi_t = phi.p(), phi2_t = phi_t * phi_t;

    return 2. * zpm / vec_qt.pt2() *
           ((osp_.term_ll * 4. * zpm * zpm * q2 * phi2_0) +
            (osp_.term_tt1 * (zp * zp + zm * zm) * phi2_t + mf2_ * phi2_0) +
            (osp_.term_tt2 * (cross * cross - dot * dot) / vec_pho.pt2()) -
            (osp_.term_lt * 4. * zpm * (zp - zm) * phi_0 * dot));
  };
  //--- t-channel
  const auto q2_1 = utils::kt::q2(x1(), q1().pt2(), mA2(), mX2());
  const auto [zp_1, zm_1] = compute_zs(+1, x1());
  const auto amat2_1 = compute_mat_element(zp_1, zm_1, q2_1, q1(), q2().transverse());

  //--- u-channel
  const auto q2_2 = utils::kt::q2(x2(), q2().pt2(), mB2(), mY2());
  const auto [zp_2, zm_2] = compute_zs(-1, x2());
  const auto amat2_2 = compute_mat_element(zp_2, zm_2, q2_2, q2(), q1().transverse());

  //--- symmetrisation
  const auto amat2 = 0.5 * (osp_.mat1 * amat2_1 + osp_.mat2 * amat2_2);
  if (!utils::positive(amat2))
    return 0.;

  const auto t_limits = Limits{0., std::pow(std::max(mt1, mt2), 2)};
  const auto prefac = g_part1_(std::sqrt(t_limits.trim(q2_1))) * g_part2_(std::sqrt(t_limits.trim(q2_2)));
  if (!utils::positive(prefac))
    return 0.;
  return prefac * amat2;
}
// register process
REGISTER_PROCESS("pptoff", PPtoFF);
