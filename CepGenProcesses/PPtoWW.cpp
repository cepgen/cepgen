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
#include "CepGen/Physics/PDG.h"
#include "CepGen/Processes/Process2to4.h"

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
      double offShellME(double phi_sum, double phi_diff) const;

      double amplitudeWW(double shat, double that, double uhat, short lam1, short lam2, short lam3, short lam4) const;
      enum class AmplitudesMode { SM, W, Wbar, phiW, WB };

      static constexpr double prefactor_ = constants::G_EM_SQ * constants::G_EM_SQ;

      const double mW_, mW2_;
      const int method_;
      struct EftParameters {
        explicit EftParameters(const ParametersList& params)
            : mode(params.getAs<int, AmplitudesMode>("mode", AmplitudesMode::SM)),
              s1(params.get<double>("s1")),
              mH(params.get<double>("mH")) {}
        AmplitudesMode mode;
        double s1, mH;
        double c1() const { return sqrt(1. - s1 * s1); }
      } eft_ext_;

      std::vector<short> pol_w1_, pol_w2_;
    };

    PPtoWW::PPtoWW(const ParametersList& params)
        : Process2to4(params, {PDG::photon, PDG::photon}, PDG::W),
          mW_(PDG::get().mass(PDG::W)),
          mW2_(mW_ * mW_),
          method_(params.get<int>("method", 1)),
          eft_ext_(params.get<ParametersList>("eftExtension")) {
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
      CG_DEBUG("PPtoWW:mode") << "matrix element computation method: " << method_ << ".";

      if (eft_ext_.mode != AmplitudesMode::SM)
        CG_INFO("PPtoWW") << "EFT extension enabled. Parameters: " << params.get<ParametersList>("eftExtension") << ".";
    }

    void PPtoWW::prepareProcessKinematics() {
      cuts::Central single_w_cuts(cepgen::ParametersList{});
      if (kin_.cuts().central_particles.count(PDG::W) > 0)
        single_w_cuts = kin_.cuts().central_particles.at(PDG::W);
      setCuts(single_w_cuts);
    }

    double PPtoWW::computeCentralMatrixElement() const {
      CG_DEBUG_LOOP("PPtoWW:ME") << "matrix element mode: " << method_ << ".";

      double mat_el = prefactor_;
      switch (method_) {
        case 0: {
          // On-shell matrix element
          // references:
          //  Phys.Rev.D 51 (1995) 4738
          //  JHEP 02 (2015) 098
          mat_el *= onShellME();
        } break;
        case 1: {
          mat_el *= offShellME(phi_qt1_ + phi_qt2_, phi_qt1_ - phi_qt2_);
        } break;
        default:
          throw CG_FATAL("PPtoWW:ME") << "Invalid ME calculation method (" << method_ << ")!";
      }
      CG_DEBUG_LOOP("PPtoWW:ME") << "prefactor: " << prefactor_ << "\n\t"
                                 << "matrix element: " << mat_el << ".";
      return mat_el;
    }

    double PPtoWW::onShellME() const {
      const double s_hat = shat(), t_hat = that(), u_hat = uhat();

      const double term1 = 2. * s_hat * (2. * s_hat + 3. * mW2_) / (3. * (mW2_ - t_hat) * (mW2_ - u_hat));
      const double term2 =
          2. * s_hat * s_hat * (s_hat * s_hat + 3. * mW2_ * mW2_) / (3. * pow(mW2_ - t_hat, 2) * pow(mW2_ - u_hat, 2));

      return 6. * (1. - term1 + term2);
    }

    double PPtoWW::offShellME(double phi_sum, double phi_diff) const {
      const double s_hat = shat(), t_hat = that(), u_hat = uhat();
      double amat2_0 = 0., amat2_1 = 0., amat2_interf = 0.;
      for (const auto lam3 : pol_w1_)
        for (const auto lam4 : pol_w2_) {
          double ampli_pp = amplitudeWW(s_hat, t_hat, u_hat, +1, +1, lam3, lam4);
          double ampli_mm = amplitudeWW(s_hat, t_hat, u_hat, -1, -1, lam3, lam4);
          double ampli_pm = amplitudeWW(s_hat, t_hat, u_hat, +1, -1, lam3, lam4);
          double ampli_mp = amplitudeWW(s_hat, t_hat, u_hat, -1, +1, lam3, lam4);

          amat2_0 += ampli_pp * ampli_pp + ampli_mm * ampli_mm + 2. * cos(2. * phi_diff) * ampli_pp * ampli_mm;
          amat2_1 += ampli_pm * ampli_pm + ampli_mp * ampli_mp + 2. * cos(2. * phi_sum) * ampli_pm * ampli_mp;
          amat2_interf -= 2. * (cos(phi_sum + phi_diff) * (ampli_pp * ampli_pm + ampli_mm * ampli_mp) +
                                cos(phi_sum - phi_diff) * (ampli_pp * ampli_mp + ampli_mm * ampli_pm));
        }
      return amat2_0 + amat2_1 + amat2_interf;
    }

    double PPtoWW::amplitudeWW(
        double shat, double that, double uhat, short lam1, short lam2, short lam3, short lam4) const {
      //--- first compute some kinematic variables
      const double beta2 = 1. - 4. * mW2_ / shat, beta = sqrt(beta2);
      const double inv_gamma2 = 1. - beta2, gamma2 = 1. / inv_gamma2;
      const double gamma = sqrt(gamma2), inv_gamma = 1. / gamma;
      const double cos_theta = (that - uhat) / shat / beta, cos_theta2 = cos_theta * cos_theta;
      const double sin_theta2 = 1. - cos_theta2, sin_theta = sqrt(sin_theta2);
      const double invA = 1. / (1. - beta2 * cos_theta2);

      //--- per-helicity amplitude

      switch (eft_ext_.mode) {
        case AmplitudesMode::SM: {
          if (lam3 == 0 && lam4 == 0)  // longitudinal-longitudinal
            return invA * inv_gamma2 * ((gamma2 + 1.) * (1. - lam1 * lam2) * sin_theta2 - (1. + lam1 * lam2));

          if (lam4 == 0)  // transverse-longitudinal
            return invA * (-M_SQRT2 * inv_gamma * (lam1 - lam2) * (1. + lam1 * lam3 * cos_theta) * sin_theta);

          if (lam3 == 0)  // longitudinal-transverse
            return invA * (-M_SQRT2 * inv_gamma * (lam2 - lam1) * (1. + lam2 * lam4 * cos_theta) * sin_theta);

          // transverse-transverse
          return -0.5 * invA *
                 (2. * beta * (lam1 + lam2) * (lam3 + lam4) -
                  inv_gamma2 * (1. + lam3 * lam4) * (2. * lam1 * lam2 + (1. - lam1 * lam2) * cos_theta2) +
                  (1. + lam1 * lam2 * lam3 * lam4) * (3. + lam1 * lam2) +
                  2. * (lam1 - lam2) * (lam3 - lam4) * cos_theta +
                  (1. - lam1 * lam2) * (1. - lam3 * lam4) * cos_theta2);
        }
        case AmplitudesMode::W: {
          if (lam3 == 0 && lam4 == 0)  // longitudinal-longitudinal
            return 3. * shat * eft_ext_.s1 * M_SQRT2 * constants::G_F * invA * inv_gamma2 * sin_theta2 *
                   (1. + lam1 * lam2);

          if (lam4 == 0)  // transverse-longitudinal
            return 1.5 * shat * eft_ext_.s1 * constants::G_F * invA * inv_gamma * sin_theta *
                   ((lam1 - lam2) * beta2 - beta * cos_theta * (lam1 + lam2) -
                    2 * lam3 * cos_theta * (lam1 * lam2 + inv_gamma2));

          if (lam3 == 0)  // longitudinal-transverse
            return 1.5 * shat * eft_ext_.s1 * constants::G_F * invA * inv_gamma * sin_theta *
                   ((lam2 - lam1) * beta2 - beta * cos_theta * (lam2 + lam1) -
                    2 * lam4 * cos_theta * (lam2 * lam1 + inv_gamma2));

          // transverse-transverse
          return 0.75 * shat * eft_ext_.s1 * M_SQRT2 * constants::G_F *
                 (-inv_gamma2 * beta * (1 + cos_theta2) * (lam1 + lam2) * (lam3 + lam4) +
                  2 * sin_theta2 *
                      (3. + lam3 * lam4 + lam1 * lam2 * (1 - lam3 * lam4) - beta * (lam1 + lam2) * (lam3 + lam4)) -
                  2 * inv_gamma2 *
                      (2 + (1 - lam1 * lam2) * lam3 * lam4 - cos_theta2 * (3 + lam1 * lam2 + 2 * lam3 * lam4)));
        }
        case AmplitudesMode::Wbar: {
          if (lam3 == 0 && lam4 == 0)  // longitudinal-longitudinal
            return -3 * shat * eft_ext_.s1 * M_SQRT2 * constants::G_F * inv_gamma2 * invA * sin_theta2 * (lam1 + lam2);

          if (lam4 == 0)  // transverse-longitudinal
            return 1.5 * shat * eft_ext_.s1 * constants::G_F * inv_gamma * invA * sin_theta *
                   (beta * (lam1 - lam2) * lam3 + cos_theta * (2 * beta + (2. - beta2) * (lam1 + lam2) * lam3));

          if (lam3 == 0)  // longitudinal-transverse
            return 1.5 * shat * eft_ext_.s1 * constants::G_F * inv_gamma * invA * sin_theta *
                   (beta * (lam2 - lam1) * lam4 + cos_theta * (2 * beta + (2. - beta2) * (lam2 + lam1) * lam4));

          // transverse-transverse
          return -1.5 * shat * eft_ext_.s1 * M_SQRT2 * constants::G_F * invA *
                 (2 * sin_theta2 * (lam1 + lam2 - beta * (lam3 + lam4)) +
                  inv_gamma2 * ((lam1 + lam2) * (cos_theta2 * (2 + lam3 * lam4) - 1) -
                                beta * (cos_theta2 + lam1 * lam2) * (lam3 + lam4)));
        }
        case AmplitudesMode::phiW: {
          const double invB = 1. / (shat - eft_ext_.mH * eft_ext_.mH);
          if (lam3 == 0 && lam4 == 0)  // longitudinal-longitudinal
            return -0.25 * shat * shat * eft_ext_.s1 * eft_ext_.s1 * M_SQRT2 * constants::G_F * invB * (1 + beta2) *
                   (1 + lam1 * lam2);

          if (lam4 == 0 || lam3 == 0)  // transverse-longitudinal or longitudinal-transverse
            return 0.;

          // transverse-transverse
          return -0.125 * shat * shat * eft_ext_.s1 * eft_ext_.s1 * M_SQRT2 * constants::G_F * inv_gamma2 * invB *
                 (1 + lam1 * lam2) * (1 + lam3 * lam4);
        }
        case AmplitudesMode::WB: {
          const double invB = 1. / (shat - eft_ext_.mH * eft_ext_.mH);
          if (lam3 == 0 && lam4 == 0)  // longitudinal-longitudinal
            return 2 * invA * eft_ext_.c1() / eft_ext_.s1 *
                       (1 - lam1 * lam2 - 2 * cos_theta2 - gamma2 * (1 + lam1 * lam2) * sin_theta2) +
                   0.5 * shat * shat * constants::G_F * M_SQRT2 * invB * eft_ext_.s1 * eft_ext_.c1() * (1 + beta2) *
                       (1 + lam1 * lam2);

          if (lam4 == 0)  // transverse-longitudinal
            return 0.5 * gamma * M_SQRT2 * invB * eft_ext_.c1() / eft_ext_.s1 * sin_theta *
                   ((lam2 - lam1) * (1 + inv_gamma2) +
                    (beta * (lam1 + lam2) + 2 * lam3 * (lam1 * lam2 - inv_gamma2)) * cos_theta);

          if (lam3 == 0)  // longitudinal-transverse
            return 0.5 * gamma * M_SQRT2 * invB * eft_ext_.c1() / eft_ext_.s1 * sin_theta *
                   ((lam1 - lam2) * (1 + inv_gamma2) +
                    (beta * (lam2 + lam1) + 2 * lam4 * (lam2 * lam1 - inv_gamma2)) * cos_theta);

          // transverse-transverse
          return -0.5 * invA * eft_ext_.c1() / eft_ext_.s1 *
                     (beta * (lam1 + lam2) * (lam3 + lam4) * (1 + cos_theta2) +
                      2 * (2 + (lam1 - lam2) * (lam3 - lam4) * cos_theta +
                           ((lam1 * lam2 - 1) * cos_theta2 + 1 + lam1 * lam2) * lam3 * lam4)) +
                 0.25 * shat * shat * M_SQRT2 * constants::G_F * inv_gamma2 * invB * eft_ext_.s1 * eft_ext_.c1() *
                     (1 + lam1 * lam2) * (1 + lam3 * lam4);
        }
      }
      throw CG_FATAL("PPtoWW:WWamplitudes") << "Invalid mode: " << (int)eft_ext_.mode << "!";
    }
  }  // namespace proc
}  // namespace cepgen
// register process
REGISTER_PROCESS("pptoww", PPtoWW)
