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
#include "CepGen/HelicityAmplitudes/PolarisationState.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "PPtoPair.h"

using namespace cepgen;

/// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow ZZ\f$ process using \f$k_{\rm T}\f$-factorization approach
class PPtoZZ final : public proc::PPtoPair {
public:
  explicit PPtoZZ(const ParametersList& params) : PPtoPair(params, PolarisationState{polarisationParameters(params)}) {}

  proc::ProcessPtr clone() const override { return std::make_unique<PPtoZZ>(parameters()); }

  static ParametersDescription description() {
    auto desc = PPtoPair::description();
    desc.setDescription("γγ → ZZ (anomalous couplings only)");
    desc.addAs<int, pdgid_t>("pair", PDG::Z).setDescription("type of central particles emitted");
    desc.add("ktFactorised", true);
    desc.addAs<int, Mode>("method", Mode::helicity)
        .setDescription("Matrix element computation method")
        .allow(2, "off-shell by Nachtmann et al. (helicity amplitudes)");
    desc.add("model", amplitudes::HelicityAmplitudesFactory::get().describeParameters("dim6EFTzz"))
        .setDescription("Helicity amplitudes model");
    desc.add("polarisationStates", PolarisationState::description());
    return desc;
  }

private:
  /// Build the polarisation states steering, also accepting the integer mode used by steering cards
  /// \note Static, as it is evaluated before the base classes are constructed.
  static ParametersList polarisationParameters(const ParametersList& params) {
    if (params.has<int>("polarisationStates"))
      return ParametersList().set<int>("mode", params.get<int>("polarisationStates"));
    return params.get<ParametersList>("polarisationStates");
  }

  void prepareFactorisedPhaseSpace() override {
    if (method_ == Mode::offShell) {
      CG_WARNING("PPtoZZ:ME") << "Method 1 (off-shell by Nachtmann et al.) is now referred to as method 2. "
                              << "Please update your steering.";
      method_ = Mode::helicity;
    }
    PPtoPair::prepareFactorisedPhaseSpace();
    if (method_ == Mode::helicity) {  // invalid steering fails at initialisation
      if (!amplitudes_->hasAAZZ())
        throw CG_FATAL("PPtoZZ:prepare") << "Invalid EFT extension enabled for γγ → ZZ (AAZZ). Model '"
                                         << amplitudes_->name()
                                         << "' does not have a registered AAZZ vertex in its specifications.";
      CG_INFO("PPtoZZ:prepare") << "Helicity amplitudes (model: " << amplitudes_->name()
                                << ") initialised. Parameters: " << amplitudes_->parameters() << ".";
    }
    CG_DEBUG("PPtoZZ:prepare") << "matrix element computation method: " << static_cast<int>(method_) << ", "
                               << "polarisation states: W1=" << polarisation_state_.polarisations().first
                               << ", W2=" << polarisation_state_.polarisations().second << ".";
  }
  double computeFactorisedMatrixElement() override {
    CG_DEBUG_LOOP("PPtoZZ:ME") << "matrix element mode: " << static_cast<int>(method_) << ".";
    switch (method_) {
      case Mode::helicity:
        return helicityME();
      default:
        throw CG_FATAL("PPtoZZ:ME") << "Invalid ME calculation method (" << static_cast<int>(method_) << ")!";
    }
  }
  double helicityPrefactor(double, double) const override { return 1.; }  ///< W+W- ampltiudes already have em prefactor
};
REGISTER_PROCESS("pptozz", PPtoZZ);
