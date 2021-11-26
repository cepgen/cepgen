/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include "CepGen/Processes/Process2to4.h"

using namespace std::complex_literals;

namespace cepgen {
  namespace proc {
    /// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process using \f$k_{\rm T}\f$-factorization approach
    /// \note The full theoretical description of this process definition may be found in \cite Luszczak:2018ntp.
    class PPtoWW final : public Process2to4 {
    public:
      explicit PPtoWW(const ParametersList&);
      ProcessPtr clone() const override { return ProcessPtr(new PPtoWW(*this)); }
      enum class Polarisation { full = 0, LL = 1, LT = 2, TL = 3, TT = 4 };
      static std::string description() { return "ɣɣ → W⁺W¯ (kt-factor.)"; }

    private:
      void prepareProcessKinematics() override;
      double computeCentralMatrixElement() const override;

      double onShellME() const;
      double offShellME() const;

      const double mW_, mW2_;
      const int method_;
      NachtmannAmplitudes ampl_;

      std::vector<int> pol_w1_, pol_w2_;
    };

    PPtoWW::PPtoWW(const ParametersList& params)
        : Process2to4(params, {PDG::photon, PDG::photon}, PDG::W),
          mW_(PDG::get().mass(PDG::W)),
          mW2_(mW_ * mW_),
          method_(params.get<int>("method", 1)),
          ampl_(params) {
      if (params.has<int>("polarisationStates"))
        switch (params.getAs<int, Polarisation>("polarisationStates", Polarisation::full)) {
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
      else if (params.has<ParametersList>("polarisationStates")) {
        const auto& states = params.get<ParametersList>("polarisationStates");
        pol_w1_ = states.get<std::vector<int> >("W1");
        pol_w2_ = states.get<std::vector<int> >("W2");
      }
      CG_DEBUG("PPtoWW") << "matrix element computation method: " << method_ << ", "
                         << "polarisation states: W1=" << pol_w1_ << ", W2=" << pol_w2_ << ".";

      if (method_ == 1) {
        CG_INFO("PPtoWW") << "Nachtmann amplitudes (model: " << ampl_.mode() << ") initialised.";
        if (ampl_.mode() != NachtmannAmplitudes::Mode::SM) {
          if (ampl_.mode() != NachtmannAmplitudes::Mode::W && ampl_.mode() != NachtmannAmplitudes::Mode::Wbar)
            throw CG_FATAL("PPtoWW") << "Invalid EFT extension enabled for ɣɣ → W⁺W¯! "
                                     << "Only supported extensions are W and Wbar. Specified model: " << ampl_.mode()
                                     << ".";
          CG_INFO("PPtoWW") << "EFT extension enabled. Parameters: " << params.get<ParametersList>("eftParameters")
                            << ".";
        }
      }
    }

    void PPtoWW::prepareProcessKinematics() {
      cuts::Central single_w_cuts(cepgen::ParametersList{});
      if (kin_.cuts().central_particles.count(PDG::W) > 0)
        single_w_cuts = kin_.cuts().central_particles.at(PDG::W);
      setCuts(single_w_cuts);
    }

    double PPtoWW::computeCentralMatrixElement() const {
      CG_DEBUG_LOOP("PPtoWW:ME") << "matrix element mode: " << method_ << ".";

      double mat_el = 1.;
      switch (method_) {
        case 0: {
          // On-shell matrix element
          // references:
          //  Phys.Rev.D 51 (1995) 4738
          //  JHEP 02 (2015) 098
          mat_el = onShellME();
        } break;
        case 1: {
          mat_el = onShellME();
        } break;
        default:
          throw CG_FATAL("PPtoWW:ME") << "Invalid ME calculation method (" << method_ << ")!";
      }
      CG_DEBUG_LOOP("PPtoWW:ME") << "matrix element: " << mat_el << ".";
      return mat_el;
    }

    double PPtoWW::onShellME() const {
      const double s_hat = shat(), t_hat = that(), u_hat = uhat();

      const double term1 = 2. * s_hat * (2. * s_hat + 3. * mW2_) / (3. * (mW2_ - t_hat) * (mW2_ - u_hat));
      const double term2 =
          2. * s_hat * s_hat * (s_hat * s_hat + 3. * mW2_ * mW2_) / (3. * pow(mW2_ - t_hat, 2) * pow(mW2_ - u_hat, 2));

      return 6. * constants::G_EM_SQ * constants::G_EM_SQ * (1. - term1 + term2);
    }

    double PPtoWW::offShellME() const {
      const NachtmannAmplitudes::Kinematics kin(mW2_, shat(), that(), uhat());
      double hel_mat_elem{0.};
      for (const auto& lam3 : pol_w1_)
        for (const auto& lam4 : pol_w2_) {
          const auto ampli_pp = ampl_(kin, +1, +1, lam3, lam4);
          const auto ampli_mm = ampl_(kin, -1, -1, lam3, lam4);
          const auto ampli_pm = ampl_(kin, +1, -1, lam3, lam4);
          const auto ampli_mp = ampl_(kin, -1, +1, lam3, lam4);

          hel_mat_elem += norm(((q1_.px() * q2_.px() + q1_.py() * q2_.py()) * (ampli_pp + ampli_mm) -
                                1.i * (q1_.px() * q2_.py() - q1_.py() * q2_.px()) * (ampli_pp - ampli_mm) -
                                (q1_.px() * q2_.px() - q1_.py() * q2_.py()) * (ampli_pm + ampli_mp) -
                                1.i * (q1_.px() * q2_.py() + q1_.py() * q2_.px()) * (ampli_pm - ampli_mp)) *
                               0.5 / qt1_ / qt2_);
        }
      return hel_mat_elem;
    }
  }  // namespace proc
}  // namespace cepgen
// register process
REGISTER_PROCESS("pptoww", PPtoWW)
