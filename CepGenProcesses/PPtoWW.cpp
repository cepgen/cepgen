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

#include <complex>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"
#include "CepGen/HelicityAmplitudes/PolarisationState.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGenProcesses/PPtoVV.h"

using namespace cepgen;
using namespace std::string_literals;

/// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process using \f$k_{\rm T}\f$-factorization approach
/// \note The full theoretical description of this process definition may be found in \cite Luszczak:2018ntp.
class PPtoWW final : public PPtoVV {
public:
  explicit PPtoWW(const ParametersList& params)
      : PPtoVV(params, {+static_cast<spdgid_t>(PDG::W), -static_cast<spdgid_t>(PDG::W)}, PDG::W),
        method_(steer<int>("method")) {
    CG_DEBUG("PPtoWW") << "matrix element computation method: " << method_ << ", "
                       << "polarisation states: "
                       << "W1=" << polarisation_.polarisations().first << ", "
                       << "W2=" << polarisation_.polarisations().second << ".";
    if (method_ == 1) {  // build the amplitudes model once, so that invalid steering fails at initialisation
      const auto& model_params = steer<ParametersList>("model");
      if (amplitudes_ = amplitudes::HelicityAmplitudesFactory::get().build(model_params); !amplitudes_->hasAAWW())
        throw CG_FATAL("PPtoWW") << "Invalid EFT extension enabled for γγ → W⁺W¯ (AAWW). Model '" << amplitudes_->name()
                                 << "' does not have a registered AAWW vertex in its specifications.";
      CG_INFO("PPtoWW") << "Helicity amplitudes (model: " << amplitudes_->name()
                        << ") initialised. Parameters: " << model_params << ".";
    }
  }

  proc::ProcessPtr clone() const override { return std::make_unique<PPtoWW>(parameters()); }

  static ParametersDescription description() {
    auto desc = PPtoVV::description();
    desc.setDescription("γγ → W⁺W¯");
    desc.add("method", 1)
        .setDescription("Matrix element computation method")
        .allow(0, "on-shell")
        .allow(1, "off-shell by Nachtmann et al.");
    desc.add("model", amplitudes::HelicityAmplitudesFactory::get().describeParameters("standardModel"));
    return desc;
  }

private:
  void prepareFactorisedPhaseSpace() override {
    if (method_ == 1 && !partonsKtFactorised())
      throw CG_FATAL("PPtoWW") << "Off-shell matrix element only defined for processes with initial partons kT.";
  }
  double computeFactorisedMatrixElement() override {
    CG_DEBUG_LOOP("PPtoWW:ME") << "matrix element mode: " << method_ << ".";
    switch (method_) {
      case 0:
        return onShellME();
      case 1:
        return offShellME();
      default:
        throw CG_FATAL("PPtoWW:ME") << "Invalid ME calculation method (" << method_ << ")!";
    }
  }

  double onShellME() const {
    // On-shell matrix element
    // references:
    //  Phys.Rev.D 51 (1995) 4738
    //  JHEP 02 (2015) 098
    const double s_hat = shat(), t_hat = that(), u_hat = uhat();

    const double term1 = 2. * s_hat * (2. * s_hat + 3. * mv2_) / (3. * (mv2_ - t_hat) * (mv2_ - u_hat));
    const double term2 = 2. * s_hat * s_hat * (s_hat * s_hat + 3. * mv2_ * mv2_) /
                         (3. * std::pow(mv2_ - t_hat, 2) * std::pow(mv2_ - u_hat, 2));

    return 6. * constants::G_EM_SQ * constants::G_EM_SQ * (1. - term1 + term2) / s_hat / s_hat;
  }

  const int method_;
};
REGISTER_PROCESS("pptoww", PPtoWW);
