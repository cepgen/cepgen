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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/NachtmannAmplitudes.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  NachtmannAmplitudes::NachtmannAmplitudes(const ParametersList& params)
      : mode_(params.getAs<int, NachtmannAmplitudes::Mode>("model", NachtmannAmplitudes::Mode::SM)),
        eft_ext_(params.get<ParametersList>("eftParameters")),
        mw2_(PDG::get().mass(24)) {
    CG_DEBUG("NachtmannAmplitudes") << "Nachtmann amplitudes evaluation framework built for mode=" << (int)mode_ << ".";
  }

  NachtmannAmplitudes::EFTParameters::EFTParameters(const ParametersList& params)
      : s1(params.get<double>("s1")), mH(params.get<double>("mH")) {}

  double NachtmannAmplitudes::operator()(
      double shat, double that, double uhat, short lam1, short lam2, short lam3, short lam4) const {
    //--- first compute some kinematic variables
    const double beta2 = 1. - 4. * mw2_ / shat, beta = sqrt(beta2);
    const double inv_gamma2 = 1. - beta2, gamma2 = 1. / inv_gamma2;
    const double gamma = sqrt(gamma2), inv_gamma = 1. / gamma;
    const double cos_theta = (that - uhat) / shat / beta, cos_theta2 = cos_theta * cos_theta;
    const double sin_theta2 = 1. - cos_theta2, sin_theta = sqrt(sin_theta2);
    const double invA = 1. / (1. - beta2 * cos_theta2);

    //--- per-helicity amplitude

    switch (mode_) {
      case Mode::SM: {
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
                (1. + lam1 * lam2 * lam3 * lam4) * (3. + lam1 * lam2) + 2. * (lam1 - lam2) * (lam3 - lam4) * cos_theta +
                (1. - lam1 * lam2) * (1. - lam3 * lam4) * cos_theta2);
      }
      case Mode::W: {
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
      case Mode::Wbar: {
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
      case Mode::phiW: {
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
      case Mode::WB: {
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
    throw CG_FATAL("PPtoWW:WWamplitudes") << "Invalid mode: " << mode_ << "!";
  }

  std::ostream& operator<<(std::ostream& os, const NachtmannAmplitudes::Mode& mode) {
    switch (mode) {
      case NachtmannAmplitudes::Mode::SM:
        return os << "Standard model";
      case NachtmannAmplitudes::Mode::W:
        return os << "W";
      case NachtmannAmplitudes::Mode::Wbar:
        return os << "W-bar";
      case NachtmannAmplitudes::Mode::phiW:
        return os << "phiW";
      case NachtmannAmplitudes::Mode::WB:
        return os << "WB";
    }
    return os << (int)mode;
  }
}  // namespace cepgen
