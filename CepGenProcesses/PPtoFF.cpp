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
#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "PPtoPair.h"

using namespace cepgen;
using namespace std::string_literals;

/// Compute the 2-to-4 matrix element for a CE \f$\gamma\gamma\rightarrow f\bar f\f$ process
class PPtoFF final : public cepgen::proc::PPtoPair {
public:
  explicit PPtoFF(const ParametersList& params) : PPtoPair(params, {-1, 1}, {-1, 1}) {}

  proc::ProcessPtr clone() const override { return std::make_unique<PPtoFF>(parameters()); }

  static ParametersDescription description() {
    auto desc = PPtoPair::description();
    desc.setDescription("γγ → f⁺f¯");
    desc.addAs<int, pdgid_t>("pair", PDG::muon).setDescription("type of central particles emitted");
    desc.add("model", amplitudes::HelicityAmplitudesFactory::get().describeParameters("fermionQED"))
        .setDescription("Helicity amplitudes model (only used for the method 2)");
    return desc;
  }

private:
  void prepareFactorisedPhaseSpace() override {
    PPtoPair::prepareFactorisedPhaseSpace();
    const auto cs_prop = PDG::get()(central().at(0));
    // define central particle properties and couplings with partons
    if (!cs_prop.fermion || cs_prop.charges.empty())
      throw CG_FATAL("PPtoFF:prepare") << "Invalid fermion pair selected: " << cs_prop << ".";
    if (method_ == Mode::helicity) {
      if (!amplitudes_->hasAAFF())
        throw CG_FATAL("PPtoFF") << "Model '" << amplitudes_->name()
                                 << "' does not have a registered γγ → f f̄ vertex in its specifications.";
      CG_INFO("PPtoFF:prepare") << "Helicity amplitudes (model: " << amplitudes_->name() << ") initialised. "
                                << "Parameters: " << amplitudes_->parameters() << ".";
    }
    CG_DEBUG("PPtoFF:prepare") << "Incoming beams: mA = " << mA() << " GeV/mB = " << mB() << " GeV.\n\t"
                               << "Produced particles: " << central() << ".\n\t"
                               << "ME computation method: " << static_cast<int>(method_) << ".";
    // constrain central particles cuts
    if (!kinematics().cuts().central.pt_diff.valid())
      kinematics().cuts().central.pt_diff = {0., 50.};
  }
  double computeFactorisedMatrixElement() override {
    switch (method_) {
      case Mode::onShell: {
        const double s_hat = shat(), t_hat = that(), u_hat = uhat();  // buffer Mandelstam variables
        CG_DEBUG_LOOP("PPtoFF:onShell") << "that: " << t_hat << ", uhat: " << u_hat << ".";

        if (t_hat == mx2_ || u_hat == mx2_)
          return 0.;
        const auto q = std::sqrt(std::fabs(t_hat));  // t-hat negative in region of interest
        const auto prefactor = colour_factor_ * g_part1_(q) * g_part2_(q);
        if (!utils::positive(prefactor))
          return 0.;

        const auto mf4 = mx2_ * mx2_, mf8 = mf4 * mf4;
        const auto out = 6. * mf8 + (-3. * mf4 * t_hat * t_hat) + (-14. * mf4 * t_hat * u_hat) +
                         (-3. * mf4 * u_hat * u_hat) + (1. * mx2_ * t_hat * t_hat * t_hat) +
                         (7. * mx2_ * t_hat * t_hat * u_hat) + (7. * mx2_ * t_hat * u_hat * u_hat) +
                         (1. * mx2_ * u_hat * u_hat * u_hat) + (-1. * t_hat * t_hat * t_hat * u_hat) +
                         (-1. * t_hat * u_hat * u_hat * u_hat);
        return -2. * prefactor * out * std::pow((mx2_ - t_hat) * (mx2_ - u_hat) * s_hat, -2);
      }
      case Mode::offShell:
        return offShellME();
      case Mode::helicity:
        return helicityME();
      default:
        throw CG_FATAL("PPtoFF") << "Invalid ME calculation method (" << static_cast<int>(method_) << ")!";
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

    return 2. * zpm / vec_qt.pt2() *
           ((osp_.term_ll * 4. * zpm * zpm * q2 * phi2_0) +
            (osp_.term_tt1 * ((zp * zp + zm * zm) * phi2_t + mx2_ * phi2_0)) +
            // azimuthal correlation term: spin-summed splitting for a linearly polarised photon gives
            // (zp^2+zm^2)|Phi|^2 + 2 zp zm (|Phi_perp|^2-|Phi_par|^2) = (zp-zm)^2 |Phi_par|^2 + |Phi_perp|^2
            (osp_.term_tt2 * 2. * zpm * (cross * cross - dot * dot) / vec_pho.pt2()) -
            (osp_.term_lt * 4. * zpm * (zp - zm) * phi_0 * dot));
  }
};

REGISTER_PROCESS("pptoff", PPtoFF);
