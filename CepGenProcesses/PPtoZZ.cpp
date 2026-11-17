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
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGenProcesses/PPtoVV.h"

using namespace cepgen;
using namespace std::string_literals;

/// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow ZZ\f$ process using \f$k_{\rm T}\f$-factorization approach
class PPtoZZ final : public PPtoVV {
public:
  explicit PPtoZZ(const ParametersList& params)
      : PPtoVV(params, {static_cast<spdgid_t>(PDG::Z), static_cast<spdgid_t>(PDG::Z)}, PDG::Z, 0.5) {
    CG_DEBUG("PPtoZZ") << "polarisation states: "
                       << "Z1=" << polarisation_.polarisations().first << ", "
                       << "Z2=" << polarisation_.polarisations().second << ".";
    const auto& model_params = steer<ParametersList>("model");
    if (amplitudes_ = amplitudes::HelicityAmplitudesFactory::get().build(model_params); !amplitudes_->hasAAZZ())
      throw CG_FATAL("PPtoZZ") << "Invalid EFT extension enabled for γγ → ZZ (AAZZ). Model '" << amplitudes_->name()
                               << "' does not have a registered AAWW vertex in its specifications.";
    CG_INFO("PPtoZZ") << "Helicity amplitudes (model: " << amplitudes_->name()
                      << ") initialised. Parameters: " << model_params << ".";
  }

  proc::ProcessPtr clone() const override { return std::make_unique<PPtoZZ>(parameters()); }

  static ParametersDescription description() {
    auto desc = PPtoVV::description();
    desc.setDescription("γγ → ZZ (anomalous couplings only)");
    // the default model has vanishing couplings, i.e. no contribution unless steered explicitly
    desc.add("model", amplitudes::HelicityAmplitudesFactory::get().describeParameters("dim6EFTzz"));
    return desc;
  }

private:
  void prepareFactorisedPhaseSpace() override {
    if (!partonsKtFactorised())
      throw CG_FATAL("PPtoZZ") << "Off-shell matrix element only defined for processes with initial partons kT.";
  }
  double computeFactorisedMatrixElement() override { return offShellME(); }
};
REGISTER_PROCESS("pptozz", PPtoZZ);
