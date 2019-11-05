/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2023  Laurent Forthomme
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
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process2to4.h"

using namespace cepgen;

/// Matrix element for the \f$\gamma\gamma\rightarrow \tilde{f}^+\tilde{f}^-/\tilde{\chi}^+\tilde{\chi}^-/H^+H^-\f$ process
class PPtoSusy : public cepgen::proc::Process2to4 {
public:
  explicit PPtoSusy(const ParametersList& params)
      : Process2to4(params, PDG::invalid),
        pair_(steer<ParticleProperties>("pair")),
        mass2_(pair_.mass * pair_.mass),
        prefactor_(constants::G_EM_SQ * constants::G_EM_SQ) {
    if (pair_.pdgid != PDG::invalid && pair_.charge == 0.)
      throw CG_FATAL("PPtoSusy:prepare") << "Invalid SUSY pair selected: " << pair_ << ")!";

    CG_DEBUG("PPtoSusy:prepare") << "Produced particles: " << pair_ << " (mass = " << pair_.mass << " GeV).";
  }

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new PPtoSusy(*this)); }
  static ParametersDescription description() {
    auto desc = Process2to4::description();
    desc.setDescription("gamma,gamma --> ~l+~l-/~chi+~chi-/H+H-");
    return desc;
  }

private:
  void prepareProcessKinematics() override {
    if (!kinematics().cuts().central.pt_diff.valid())
      kinematics().cuts().central.pt_diff = {0., 50.};
  }
  double computeCentralMatrixElement() const override {
    // NOTE: only the on-shell formula is defined for the time being

    const double s_hat = shat();  // squared two-photon mass
    if (s_hat == 0.)
      return 0.;
    const auto inv_s_hat = 1. / s_hat;
    const double mass2_norm = mass2_ * inv_s_hat;
    const double beta2 = 1. - 4. * mass2_norm;
    if (beta2 < 0.)
      return 0.;
    const auto beta = std::sqrt(beta2);  // charginos/sleptons/H+- velocity in c.m. frame
    const double log_term = std::log((1. + beta) / (1. - beta));

    if (pair_.fermion)  // charginos
      return (2. * prefactor_) *
             ((1. + 4. * mass2_norm - 8. * mass2_norm * mass2_norm) * log_term - beta * (1 + 4. * mass2_norm));
    else  // sleptons/H+-
      return prefactor_ * (beta * (1 + 4. * mass2_norm) - 4. * mass2_norm * (1. - 2. * mass2_norm) * log_term);
  }

  const ParticleProperties pair_;
  const double mass2_;
  const double prefactor_;
};

// register process
REGISTER_PROCESS("pptosusy", PPtoSusy);
