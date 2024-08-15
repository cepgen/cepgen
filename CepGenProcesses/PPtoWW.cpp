/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2024  Laurent Forthomme
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

#include <cassert>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/NachtmannAmplitudes.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/PolarisationState.h"
#include "CepGen/Process/FactorisedProcess.h"

using namespace cepgen;
using namespace std::complex_literals;

/// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process using \f$k_{\rm T}\f$-factorization approach
/// \note The full theoretical description of this process definition may be found in \cite Luszczak:2018ntp.
class PPtoWW final : public cepgen::proc::FactorisedProcess {
public:
  explicit PPtoWW(const ParametersList& params)
      : FactorisedProcess(params, {+(spdgid_t)PDG::W, -(spdgid_t)PDG::W}),
        mW_(PDG::get().mass(PDG::W)),
        mW2_(mW_ * mW_),
        method_(steer<int>("method")),
        ampl_(params_),
        pol_(steer<ParametersList>("polarisationStates")) {
    CG_DEBUG("PPtoWW") << "matrix element computation method: " << method_ << ", "
                       << "polarisation states: "
                       << "W1=" << pol_.polarisations().first << ", "
                       << "W2=" << pol_.polarisations().second << ".";

    if (method_ == 1) {
      CG_INFO("PPtoWW") << "Nachtmann amplitudes (model: " << ampl_.mode() << ") initialised.";
      if (ampl_.mode() != NachtmannAmplitudes::Mode::SM) {
        if (ampl_.mode() != NachtmannAmplitudes::Mode::W && ampl_.mode() != NachtmannAmplitudes::Mode::Wbar)
          throw CG_FATAL("PPtoWW") << "Invalid EFT extension enabled for γγ → W⁺W¯! "
                                   << "Only supported extensions are W and Wbar. Specified model: " << ampl_.mode()
                                   << ".";
        CG_INFO("PPtoWW") << "EFT extension enabled. Parameters: " << steer<ParametersList>("eftParameters") << ".";
      }
    }
  }

  proc::ProcessPtr clone() const override { return std::make_unique<PPtoWW>(*this); }

  static ParametersDescription description() {
    auto desc = FactorisedProcess::description();
    desc.setDescription("γγ → W⁺W¯");
    desc.add<bool>("ktFactorised", true);
    desc.add<int>("method", 1)
        .setDescription("Matrix element computation method")
        .allow(0, "on-shell")
        .allow(1, "off-shell by Nachtmann et al.");
    desc.add<ParametersDescription>("polarisationStates", PolarisationState::description());
    desc += NachtmannAmplitudes::description();
    return desc;
  }

private:
  void prepareFactorisedPhaseSpace() override {
    cuts::Central single_w_cuts(ParametersList{});
    if (kinematics().cuts().central_particles.count(PDG::W) > 0)
      single_w_cuts = kinematics().cuts().central_particles.at(PDG::W);
    psgen_->setCentralCuts(single_w_cuts);
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
  enum class Polarisation { full = 0, LL = 1, LT = 2, TL = 3, TT = 4 };

  double onShellME() const {
    // On-shell matrix element
    // references:
    //  Phys.Rev.D 51 (1995) 4738
    //  JHEP 02 (2015) 098
    const double s_hat = shat(), t_hat = that(), u_hat = uhat();

    const double term1 = 2. * s_hat * (2. * s_hat + 3. * mW2_) / (3. * (mW2_ - t_hat) * (mW2_ - u_hat));
    const double term2 = 2. * s_hat * s_hat * (s_hat * s_hat + 3. * mW2_ * mW2_) /
                         (3. * std::pow(mW2_ - t_hat, 2) * std::pow(mW2_ - u_hat, 2));

    return 6. * constants::G_EM_SQ * constants::G_EM_SQ * (1. - term1 + term2) / s_hat / s_hat;
  }
  double offShellME() const {
    const auto kin = NachtmannAmplitudes::Kinematics(mW2_, shat(), that(), uhat());
    const double p1 = q1().px() * q2().px() + q1().py() * q2().py(), p2 = q1().px() * q2().py() - q1().py() * q2().px(),
                 p3 = q1().px() * q2().px() - q1().py() * q2().py(), p4 = q1().px() * q2().py() + q1().py() * q2().px();

    double hel_mat_elem{0.};
    // compute ME for each W helicity
    for (const auto& lam3 : pol_.polarisations().first)
      for (const auto& lam4 : pol_.polarisations().second) {
        // compute all photon helicity amplitudes
        const auto pp = ampl_(kin, +1, +1, lam3, lam4), mm = ampl_(kin, -1, -1, lam3, lam4),
                   pm = ampl_(kin, +1, -1, lam3, lam4), mp = ampl_(kin, -1, +1, lam3, lam4);
        // add ME for this W helicity to total ME
        hel_mat_elem += std::norm(p1 * (pp + mm) - 1i * p2 * (pp - mm) - p3 * (pm + mp) - 1i * p4 * (pm - mp));
      }
    return hel_mat_elem * std::pow(0.5 / q1().pt() / q2().pt() / shat(), 2);
  }

  const double mW_, mW2_;
  const int method_;
  const NachtmannAmplitudes ampl_;
  const PolarisationState pol_;
};
// register process
REGISTER_PROCESS("pptoww", PPtoWW);
