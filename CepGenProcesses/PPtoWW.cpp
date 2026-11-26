/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2026  Laurent Forthomme
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
#include "CepGen/HelicityAmplitudes/PolarisationState.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "PPtoPair.h"

using namespace cepgen;

/// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process using \f$k_{\rm T}\f$-factorization approach
/// \note The full theoretical description of this process definition may be found in \cite Luszczak:2018ntp.
class PPtoWW final : public proc::PPtoPair {
public:
  explicit PPtoWW(const ParametersList& params) : PPtoPair(params, PolarisationState{polarisationParameters(params)}) {}

  proc::ProcessPtr clone() const override { return std::make_unique<PPtoWW>(parameters()); }

  static ParametersDescription description() {
    auto desc = PPtoPair::description();
    desc.setDescription("γγ → W⁺W¯");
    desc.addAs<int, pdgid_t>("pair", PDG::W).setDescription("type of central particles emitted");
    desc.add("ktFactorised", true);
    desc.addAs<int, Mode>("method", Mode::helicity)
        .setDescription("Matrix element computation method")
        .allow(0, "on-shell")
        .allow(2, "off-shell by Nachtmann et al. (helicity amplitudes)");
    desc.add("model", amplitudes::HelicityAmplitudesFactory::get().describeParameters("standardModel"))
        .setDescription("Helicity amplitudes model (only used for the method 2)");
    desc.add("polarisationStates", PolarisationState::description());
    return desc;
  }

private:
  void prepareFactorisedPhaseSpace() override {
    if (method_ == Mode::offShell) {
      CG_WARNING("PPtoWW:ME") << "Method 1 (off-shell by Nachtmann et al.) is now referred to as method 2. "
                              << "Please update your steering.";
      method_ = Mode::helicity;
    }
    PPtoPair::prepareFactorisedPhaseSpace();
    if (method_ == Mode::helicity) {  // invalid steering fails at initialisation
      if (!amplitudes_->hasAAWW())
        throw CG_FATAL("PPtoWW:prepare") << "Invalid EFT extension enabled for γγ → W⁺W¯ (AAWW). Model '"
                                         << amplitudes_->name()
                                         << "' does not have a registered AAWW vertex in its specifications.";
      CG_INFO("PPtoWW:prepare") << "Helicity amplitudes (model: " << amplitudes_->name()
                                << ") initialised. Parameters: " << amplitudes_->parameters() << ".";
    }
    CG_DEBUG("PPtoWW:prepare") << "matrix element computation method: " << static_cast<int>(method_) << ", "
                               << "polarisation states: W1=" << polarisation_state_.polarisations().first
                               << ", W2=" << polarisation_state_.polarisations().second << ".";
  }
  double computeFactorisedMatrixElement() override {
    CG_DEBUG_LOOP("PPtoWW:ME") << "matrix element mode: " << static_cast<int>(method_) << ".";
    switch (method_) {
      case Mode::onShell: {  // on-shell matrix element (references: Phys.Rev.D 51 (1995) 4738, JHEP 02 (2015) 098)
        const double s_hat = shat(), t_hat = that(), u_hat = uhat();

        const double term1 = 2. * s_hat * (2. * s_hat + 3. * mx2_) / (3. * (mx2_ - t_hat) * (mx2_ - u_hat));
        const double term2 = 2. * s_hat * s_hat * (s_hat * s_hat + 3. * mx2_ * mx2_) /
                             (3. * std::pow(mx2_ - t_hat, 2) * std::pow(mx2_ - u_hat, 2));

        return 6. * constants::G_EM_SQ * constants::G_EM_SQ * (1. - term1 + term2) / s_hat / s_hat;
      }
      case Mode::helicity:
        return helicityME();
      default:
        throw CG_FATAL("PPtoWW:ME") << "Invalid ME calculation method (" << static_cast<int>(method_) << ")!";
    }
  }
  double helicityPrefactor(double, double) const override { return 1.; }  ///< W+W- ampltiudes already have em prefactor
};
REGISTER_PROCESS("pptoww", PPtoWW);
