/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2023  Laurent Forthomme
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
#include "CepGen/Process/Process2to4.h"

using namespace cepgen;

/// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process using \f$k_{\rm T}\f$-factorization approach
/// \note The full theoretical description of this process definition may be found in \cite Luszczak:2018ntp.
class PPtoWW final : public cepgen::proc::Process2to4 {
public:
  explicit PPtoWW(const ParametersList& params)
      : Process2to4(params, PDG::W),
        mW_(PDG::get().mass(PDG::W)),
        mW2_(mW_ * mW_),
        method_(steer<int>("method")),
        ampl_(params_) {
    const auto& states = steer<ParametersList>("polarisationStates");
    if (!states.empty()) {
      pol_w1_ = states.get<std::vector<int> >("W1");
      pol_w2_ = states.get<std::vector<int> >("W2");
    }
    if (params_.has<int>("polarisationStates"))
      switch (steerAs<int, Polarisation>("polarisationStates")) {
        case Polarisation::LL:
          pol_w1_ = {0};
          pol_w2_ = {0};
          break;
        case Polarisation::LT:
          pol_w1_ = {0};
          pol_w2_ = {-1, 1};
          break;
        case Polarisation::TL:
          pol_w1_ = {-1, 1};
          pol_w2_ = {0};
          break;
        case Polarisation::TT:
          pol_w1_ = {-1, 1};
          pol_w2_ = {-1, 1};
          break;
        case Polarisation::full:
          pol_w1_ = {-1, 0, 1};
          pol_w2_ = {-1, 0, 1};
          break;
      }
    CG_DEBUG("PPtoWW") << "matrix element computation method: " << method_ << ", "
                       << "polarisation states: W1=" << pol_w1_ << ", W2=" << pol_w2_ << ".";

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

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new PPtoWW(*this)); }

  static ParametersDescription description() {
    auto desc = Process2to4::description();
    desc.setDescription("γγ → W⁺W¯");
    desc.add<bool>("ktFactorised", true);
    desc.add<int>("method", 1)
        .setDescription("Matrix element computation method (0 = on-shell, 1 = off-shell by Nachtmann et al.)");
    ParametersDescription pol_states;
    pol_states.add<std::vector<int> >("W1", {-1, 0, 1}).setDescription("First W+- polarisation states");
    pol_states.add<std::vector<int> >("W2", {-1, 0, 1}).setDescription("Second W+- polarisation states");
    desc.add<ParametersDescription>("polarisationStates", pol_states);
    desc += NachtmannAmplitudes::description();
    return desc;
  }

private:
  void prepareProcessKinematics() override {
    cuts::Central single_w_cuts(ParametersList{});
    if (kinematics().cuts().central_particles.count(PDG::W) > 0)
      single_w_cuts = kinematics().cuts().central_particles.at(PDG::W);
    setCuts(single_w_cuts);
  }
  double computeCentralMatrixElement() const override {
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
    const double term2 =
        2. * s_hat * s_hat * (s_hat * s_hat + 3. * mW2_ * mW2_) / (3. * pow(mW2_ - t_hat, 2) * pow(mW2_ - u_hat, 2));

    return 6. * constants::G_EM_SQ * constants::G_EM_SQ * (1. - term1 + term2);
  }
  double offShellME() const {
    const NachtmannAmplitudes::Kinematics kin(mW2_, shat(), that(), uhat());
    const double p1 = q1().px() * q2().px() + q1().py() * q2().py(), p2 = q1().px() * q2().py() - q1().py() * q2().px(),
                 p3 = q1().px() * q2().px() - q1().py() * q2().py(), p4 = q1().px() * q2().py() + q1().py() * q2().px();

    double hel_mat_elem{0.};
    // compute ME for each W helicity
    for (const auto& lam3 : pol_w1_)
      for (const auto& lam4 : pol_w2_) {
        // compute all photon helicity amplitudes
        const auto pp = ampl_(kin, +1, +1, lam3, lam4), mm = ampl_(kin, -1, -1, lam3, lam4),
                   pm = ampl_(kin, +1, -1, lam3, lam4), mp = ampl_(kin, -1, +1, lam3, lam4);
        // add ME for this W helicity to total ME
        hel_mat_elem += norm(p1 * (pp + mm) - std::complex<double>(0, 1) * p2 * (pp - mm) - p3 * (pm + mp) -
                             std::complex<double>(0, 1) * p4 * (pm - mp));
      }
    return hel_mat_elem * std::pow(0.5 / q1().p() / q2().p(), 2);
  }

  const double mW_, mW2_;
  const int method_;
  const NachtmannAmplitudes ampl_;

  std::vector<int> pol_w1_, pol_w2_;
};
// register process
REGISTER_PROCESS("pptoww", PPtoWW);
