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
#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "PPtoPair.h"

using namespace cepgen;
using namespace std::string_literals;

/// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow S^+S^-\f$ process, with \f$S\f$ a charged scalar
class PPtoSS final : public cepgen::proc::PPtoPair {
public:
  explicit PPtoSS(const ParametersList& params) : PPtoPair(params, PolarisationState{PolarisationState::Mode::LL}) {}

  proc::ProcessPtr clone() const override { return std::make_unique<PPtoSS>(parameters()); }

  static ParametersDescription description() {
    auto desc = PPtoPair::description();
    desc.setDescription("γγ → S⁺S⁻ (charged scalars pair, e.g. sleptons)");
    desc.addAs<int, pdgid_t>("pair", PDG::piPlus)
        .setDescription("type of central charged scalar emitted (must be defined in the PDG library)");
    desc.add("model", amplitudes::HelicityAmplitudesFactory::get().describeParameters("scalarQED"))
        .setDescription("Helicity amplitudes model (only used for the method 2)");
    return desc;
  }

private:
  static spdgids_t make_pdgids_pair(pdgid_t pair) {
    return spdgids_t{static_cast<spdgid_t>(pair), -static_cast<spdgid_t>(pair)};
  }

  void prepareFactorisedPhaseSpace() override {
    PPtoPair::prepareFactorisedPhaseSpace();
    const auto& cs_prop = PDG::get()(central().at(0));
    if (cs_prop.fermion)
      throw CG_FATAL("PPtoSS:prepare") << "Invalid central particle selected (fermion; use 'pptoff' instead): "
                                       << cs_prop << ".";
    if (method_ == Mode::helicity)
      CG_INFO("PPtoSS:prepare") << "Helicity amplitudes (model: " << amplitudes_->name() << ") initialised. "
                                << "Parameters: " << amplitudes_->parameters() << ".";
    CG_DEBUG("PPtoSS:prepare") << "Incoming beams: mA = " << mA() << " GeV/mB = " << mB() << " GeV.\n\t"
                               << "Produced particles: " << central() << ".";
    // constrain central particles cuts
    if (!kinematics().cuts().central.pt_diff.valid())
      kinematics().cuts().central.pt_diff = {0., 50.};
  }

  double computeFactorisedMatrixElement() override {
    switch (method_) {
      case Mode::onShell: {
        const double s_hat = shat(), t_hat = that(), u_hat = uhat();            // buffer Mandelstam variables
        const auto d = 4. * (t_hat - mx2_) * (u_hat - mx2_) / (s_hat * s_hat);  // = 1-beta^2 cos^2(theta)
        if (!utils::positive(d))
          return 0.;
        const auto q = std::sqrt(std::fabs(t_hat));  // t-hat negative in region of interest
        const auto prefactor = colour_factor_ * kFourPi * kFourPi * qx2_ * qx2_ * alphaEM(q) * alphaEM(q);
        if (!utils::positive(prefactor))
          return 0.;

        const auto one_minus_beta2 = 4. * mx2_ / s_hat, beta2_sin2 = d - one_minus_beta2;
        return 2. * prefactor * (one_minus_beta2 * one_minus_beta2 + beta2_sin2 * beta2_sin2) / (d * d) /
               (s_hat * s_hat);
      }
      case Mode::offShell:
        return offShellME();
      case Mode::helicity:
        return helicityME();
      default:
        throw CG_FATAL("PPtoSS") << "Invalid ME calculation method (" << static_cast<int>(method_) << ")!";
    }
  }

  double phiMatrixElement(
      double zp, double zm, double q2, const Momentum& vec_pho, const Momentum& vec_qt) const override {
    const auto vec_kt = Momentum(zm * pc(0) - zp * pc(1)).transverse();
    const auto phi_p = vec_kt + zp * vec_qt, phi_m = vec_kt - zm * vec_qt;
    const auto zpm = zp * zm, eps2 = mx2_ + zpm * q2;

    const auto kp = 1. / (phi_p.pt2() + eps2), km = 1. / (phi_m.pt2() + eps2);
    const auto phi = Momentum(kp * phi_p - km * phi_m).setEnergy(kp - km);
    const auto dot = phi.threeProduct(vec_pho), cross = phi.crossProduct(vec_pho, Momentum::Z);

    const auto phi_0 = phi.energy(), phi2_0 = phi_0 * phi_0, phi_t = phi.p(), phi2_t = phi_t * phi_t;

    // transverse term: 4 Phi_par^2 = 2 Phi^2 + 2 (Phi_par^2-Phi_perp^2)
    return zpm * zpm / vec_qt.pt2() *
           ((osp_.term_ll * (zp - zm) * (zp - zm) * q2 * phi2_0) + (osp_.term_tt1 * 2. * phi2_t) -
            (osp_.term_tt2 * 2. * (cross * cross - dot * dot) / vec_pho.pt2()) +
            (osp_.term_lt * 4. * (zp - zm) * phi_0 * dot));
  }
};
REGISTER_PROCESS("pptoss", PPtoSS);
