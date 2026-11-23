/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;
using namespace std::string_literals;

/// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow S^+S^-\f$ process, with \f$S\f$ a charged scalar
class PPtoSS final : public proc::FactorisedProcess {
public:
  explicit PPtoSS(const ParametersList& params)
      : proc::FactorisedProcess(params, make_pdgids_pair(params.get<ParticleProperties>("pair").pdgid)),
        method_(steerAs<int, Mode>("method")),
        osp_(steer<ParametersList>("offShellParameters")) {}

  proc::ProcessPtr clone() const override { return std::make_unique<PPtoSS>(parameters()); }

  static ParametersDescription description() {
    auto desc = FactorisedProcess::description();
    desc.setDescription("γγ → S⁺S⁻ (charged scalars pair, e.g. sleptons)");
    desc.addAs<int, pdgid_t>("pair", PDG::piPlus)
        .setDescription("type of central charged scalar emitted (must be defined in the PDG library)");
    desc.addAs<int, Mode>("method", Mode::offShell)
        .setDescription("Matrix element computation method")
        .allow(0, "on-shell")
        .allow(1, "off-shell")
        .allow(2, "off-shell (helicity amplitudes)");
    desc.add("offShellParameters", OffShellParameters::description());
    desc.add("model", amplitudes::HelicityAmplitudesFactory::get().describeParameters("scalarQED"))
        .setDescription("Helicity amplitudes model (only used for the method 2)");
    return desc;
  }

private:
  static spdgids_t make_pdgids_pair(pdgid_t pair) {
    return spdgids_t{static_cast<spdgid_t>(pair), -static_cast<spdgid_t>(pair)};
  }

  void prepareFactorisedPhaseSpace() override {
    if ((method_ == Mode::offShell || method_ == Mode::helicity) && !partonsKtFactorised())
      throw CG_FATAL("PPtoSS:prepare")
          << "Off-shell matrix element only defined for factorised process with partons kt.";
    const auto& cs_prop = PDG::get()(central().at(0));
    if (cs_prop.fermion)
      throw CG_FATAL("PPtoSS:prepare") << "Invalid central particle selected (fermion; use 'pptoff' instead): "
                                       << cs_prop << ".";
    if (cs_prop.integerCharge() == 0)
      throw CG_FATAL("PPtoSS:prepare") << "The central scalar particle must be charged: " << cs_prop << ".";
    if (cs_prop.mass <= 0.)
      throw CG_FATAL("PPtoSS:prepare") << "The central scalar particle must be massive: " << cs_prop << ".";
    for (const auto role : {Particle::Role::Parton1, Particle::Role::Parton2})
      if (const auto parton_id = event().oneWithRole(role).pdgId(); parton_id != PDG::photon)
        throw CG_FATAL("PPtoSS:prepare") << "Unsupported parton id: '" << parton_id << "'. Only photons are allowed.";
    if (method_ == Mode::helicity) {
      const auto& model_params = steer<ParametersList>("model");
      if (amplitudes_ = amplitudes::HelicityAmplitudesFactory::get().build(model_params); !amplitudes_->hasAASS())
        throw CG_FATAL("PPtoSS:prepare") << "Model '" << amplitudes_->name()
                                         << "' does not have a registered γγ → S⁺S⁻ vertex in its specifications.";
      CG_INFO("PPtoSS:prepare") << "Helicity amplitudes (model: " << amplitudes_->name()
                                << ") initialised. Parameters: " << model_params << ".";
    }

    ms2_ = cs_prop.mass * cs_prop.mass;
    qs2_ = static_cast<double>(cs_prop.integerCharge() * cs_prop.integerCharge()) / 9.;
    // sum over the final state colours (scalar leptons-like particles have no colour)
    colour_factor_ = cs_prop.colours > 0 ? cs_prop.colours : 1.;

    CG_DEBUG("PPtoSS:prepare") << "Incoming beams: mA = " << mA() << " GeV/mB = " << mB() << " GeV.\n\t"
                               << "Produced particles: " << central() << ".";

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
      case Mode::helicity:
        return helicityME();
      default:
        throw CG_FATAL("PPtoSS") << "Invalid ME calculation method (" << static_cast<int>(method_) << ")!";
    }
  }

  /// Coupling factor \f$N_c(4\pi\alpha_1Q^2)(4\pi\alpha_2Q^2)\f$ for given photon scales
  double couplingFactor(double q1, double q2) const {
    return colour_factor_ * kFourPi * kFourPi * qs2_ * qs2_ * alphaEM(q1) * alphaEM(q2);
  }

  double onShellME() const {
    const double s_hat = shat(), t_hat = that(), u_hat = uhat();            // buffer Mandelstam variables
    const auto d = 4. * (t_hat - ms2_) * (u_hat - ms2_) / (s_hat * s_hat);  // = 1-beta^2 cos^2(theta)
    if (!utils::positive(d))
      return 0.;
    const auto q = std::sqrt(std::fabs(t_hat));  // t-hat negative in region of interest
    const auto prefactor = colour_factor_ * kFourPi * kFourPi * qs2_ * qs2_ * alphaEM(q) * alphaEM(q);
    if (!utils::positive(prefactor))
      return 0.;

    const auto one_minus_beta2 = 4. * ms2_ / s_hat, beta2_sin2 = d - one_minus_beta2;
    return 2. * prefactor * (one_minus_beta2 * one_minus_beta2 + beta2_sin2 * beta2_sin2) / (d * d) / (s_hat * s_hat);
  }

  double offShellME() const {
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
      const auto zpm = zp * zm, eps2 = ms2_ + zpm * q2;

      const auto kp = 1. / (phi_p.pt2() + eps2), km = 1. / (phi_m.pt2() + eps2);
      const auto phi = Momentum(kp * phi_p - km * phi_m).setEnergy(kp - km);
      const auto dot = phi.threeProduct(vec_pho), cross = phi.crossProduct(vec_pho, Momentum::Z);

      const auto phi_0 = phi.energy(), phi2_0 = phi_0 * phi_0, phi_t = phi.p(), phi2_t = phi_t * phi_t;

      // transverse term: 4 Phi_par^2 = 2 Phi^2 + 2 (Phi_par^2-Phi_perp^2)
      return zpm * zpm / vec_qt.pt2() *
             ((osp_.term_ll * (zp - zm) * (zp - zm) * q2 * phi2_0) + (osp_.term_tt1 * 2. * phi2_t) -
              (osp_.term_tt2 * 2. * (cross * cross - dot * dot) / vec_pho.pt2()) +
              (osp_.term_lt * 4. * (zp - zm) * phi_0 * dot));
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
    if (const auto prefactor = couplingFactor(std::sqrt(t_limits.trim(q2_1)), std::sqrt(t_limits.trim(q2_2)));
        utils::positive(prefactor))
      return prefactor * amat2;
    return 0.;
  }

  /// Off-shell matrix element built from the on-shell γγ → S⁺S⁻ helicity amplitudes, with the incoming photons
  /// polarised along their respective transverse momenta (see e.g. \cite Luszczak:2018ntp for the γγ → W⁺W⁻ case)
  double helicityME() const {
    if (q1().pt2() == 0. || q2().pt2() == 0.)
      return 0.;
    const auto s_hat = shat();
    if (s_hat <= 4. * ms2_)
      return 0.;
    const amplitudes::Kinematics kin(ms2_, s_hat, that(), uhat());
    // the helicity amplitudes are defined in the scattering plane: polarisations are taken relative to the pair azimuth
    const auto phi_pair = std::atan2(pc(0).py() - pc(1).py(), pc(0).px() - pc(1).px());
    const auto e1 = std::polar(1., std::atan2(q1().py(), q1().px()) - phi_pair),
               e2 = std::polar(1., std::atan2(q2().py(), q2().px()) - phi_pair);
    const auto amat2 = std::norm(amplitudes_->amplitudeLinearPolarisations(kin, 0, 0, e1, e2)) / (s_hat * s_hat);
    if (!utils::positive(amat2))
      return 0.;

    const auto mt1 = pc(0).massT(), mt2 = pc(1).massT();  // transverse masses
    const auto q2_1 = utils::kt::q2(x1(), q1().pt2(), mA2(), mX2()),
               q2_2 = utils::kt::q2(x2(), q2().pt2(), mB2(), mY2());
    const auto t_limits = Limits{0., std::pow(std::max(mt1, mt2), 2)};
    if (const auto prefactor = couplingFactor(std::sqrt(t_limits.trim(q2_1)), std::sqrt(t_limits.trim(q2_2)));
        utils::positive(prefactor))
      return prefactor * amat2;
    return 0.;
  }

  const enum class Mode { onShell = 0, offShell = 1, helicity = 2 } method_;
  std::unique_ptr<const amplitudes::HelicityAmplitudes> amplitudes_;  ///< helicity amplitudes (method 2 only)

  /// Parameters for the off-shell matrix element
  struct OffShellParameters : SteeredObject<OffShellParameters> {
    explicit OffShellParameters(const ParametersList& params) : SteeredObject(params) {
      (*this)
          .add("mat1"s, mat1)
          .add("mat2"s, mat2)
          .add("termLL"s, term_ll)
          .add("termLT"s, term_lt)
          .add("termTT"s, term_tt1)
          .add("termtt"s, term_tt2);
    }

    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.add("mat1"s, 1).setDescription("symmetry factor for the first incoming photon");
      desc.add("mat2"s, 1).setDescription("symmetry factor for the second incoming photon");
      desc.add("termLL"s, 1).setDescription("fully longitudinal relative weight");
      desc.add("termLT"s, 1).setDescription("cross-polarisation relative weight");
      desc.add("termTT"s, 1).setDescription("fully transverse (polarisation-averaged) relative weight");
      desc.add("termtt"s, 1).setDescription("transverse-transverse azimuthal correlation relative weight");
      return desc;
    }

    int mat1{0}, mat2{0};
    int term_ll{0}, term_lt{0}, term_tt1{0}, term_tt2{0};
  } osp_;

  static constexpr double kFourPi = 4. * M_PI;
  double ms2_{0.}, qs2_{0.}, colour_factor_{1.};
};
REGISTER_PROCESS("pptoss", PPtoSS);
