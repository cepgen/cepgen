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

using namespace std::complex_literals;

namespace cepgen {
  NachtmannAmplitudes::NachtmannAmplitudes(const ParametersList& params)
      : mode_(params.getAs<int, NachtmannAmplitudes::Mode>("model", NachtmannAmplitudes::Mode::SM)),
        eft_ext_(params.get<ParametersList>("eftParameters")) {
    CG_DEBUG("NachtmannAmplitudes") << "Nachtmann amplitudes evaluation framework built for mode=" << (int)mode_ << ".";
  }

  NachtmannAmplitudes::EFTParameters::EFTParameters(const ParametersList& params)
      : s1(params.get<double>("s1")), mH(params.get<double>("mH")) {}

  NachtmannAmplitudes::Kinematics::Kinematics(double mw2, double shat, double that, double uhat)
      : shat(shat),
        that(that),
        uhat(uhat),
        mw2_(mw2),
        beta2(1. - 4. * mw2_ / shat),
        beta(sqrt(beta2)),
        inv_gamma2(1. - beta2),
        gamma2(1. / inv_gamma2),
        gamma(sqrt(gamma2)),
        inv_gamma(1. / gamma),
        cos_theta((that - uhat) / shat / beta),
        cos_theta2(cos_theta * cos_theta),
        sin_theta2(1. - cos_theta2),
        sin_theta(sqrt(sin_theta2)),
        invA(1. / (1. - beta2 * cos_theta2)) {}

  bool NachtmannAmplitudes::Kinematics::operator!=(const Kinematics& oth) const {
    // checking only the base variables as all others are computed from these three
    return shat != oth.shat || that != oth.that || uhat != oth.uhat;
  }

  std::complex<double> NachtmannAmplitudes::operator()(
      const Kinematics& kin, short lam1, short lam2, short lam3, short lam4) const {
    const Helicities hel{lam1, lam2, lam3, lam4};

    //--- per-helicity amplitude
    switch (mode_) {
      case Mode::SM:
        return amplitudeSM(kin, hel);
      case Mode::W:
        return amplitudeW(kin, hel);
      case Mode::Wbar:
        return amplitudeWbar(kin, hel);
      case Mode::phiW:
        return amplitudephiW(kin, hel);
      case Mode::phiWbar:
        return 2.i * double(lam1) * amplitudephiW(kin, hel);
      case Mode::phiB:
        return std::pow(eft_ext_.c1() / eft_ext_.s1, 2) * amplitudephiW(kin, hel);
      case Mode::phiBbar:
        return 2.i * double(lam1) * std::pow(eft_ext_.c1() / eft_ext_.s1, 2) * amplitudephiW(kin, hel);
      case Mode::WB:
        return amplitudeWB(kin, hel);
      case Mode::WbarB:
        return amplitudeWbarB(kin, hel);
    }
    throw CG_FATAL("PPtoWW:WWamplitudes") << "Invalid mode: " << mode_ << "!";
  }

  std::complex<double> NachtmannAmplitudes::amplitudeSM(const Kinematics& kin, const Helicities& hel) const {
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return 1i * G_EM_SQ * kin.invA * kin.inv_gamma2 *
             ((kin.gamma2 + 1.) * (1. - hel.lam1 * hel.lam2) * kin.sin_theta2 - (1. + hel.lam1 * hel.lam2));

    if (hel.lam4 == 0)  // transverse-longitudinal
      return -1i * G_EM_SQ * kin.invA * M_SQRT2 * kin.inv_gamma * double(hel.lam1 - hel.lam2) *
             (1. + hel.lam1 * hel.lam3 * kin.cos_theta) * kin.sin_theta;

    if (hel.lam3 == 0)  // longitudinal-transverse
      return -1i * G_EM_SQ * kin.invA * M_SQRT2 * kin.inv_gamma * double(hel.lam2 - hel.lam1) *
             (1. + hel.lam2 * hel.lam4 * kin.cos_theta) * kin.sin_theta;

    // transverse-transverse
    return -0.5i * G_EM_SQ * kin.invA *
           (2. * kin.beta * float(hel.lam1 + hel.lam2) * (hel.lam3 + hel.lam4) -
            kin.inv_gamma2 * (1. + hel.lam3 * hel.lam4) *
                (2. * hel.lam1 * hel.lam2 + (1. - hel.lam1 * hel.lam2) * kin.cos_theta2) +
            (1. + hel.lam1 * hel.lam2 * hel.lam3 * hel.lam4) * (3. + hel.lam1 * hel.lam2) +
            2. * (hel.lam1 - hel.lam2) * (hel.lam3 - hel.lam4) * kin.cos_theta +
            (1. - hel.lam1 * hel.lam2) * (1. - hel.lam3 * hel.lam4) * kin.cos_theta2);
  }

  std::complex<double> NachtmannAmplitudes::amplitudeW(const Kinematics& kin, const Helicities& hel) const {
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return 3.i * G_EM * kin.shat * eft_ext_.s1 * M_SQRT2 * constants::G_F * kin.invA * kin.inv_gamma2 *
             kin.sin_theta2 * (1. + hel.lam1 * hel.lam2);

    if (hel.lam4 == 0)  // transverse-longitudinal
      return 1.5i * G_EM * kin.shat * eft_ext_.s1 * constants::G_F * kin.invA * kin.inv_gamma * kin.sin_theta *
             ((hel.lam1 - hel.lam2) * kin.beta2 - kin.beta * kin.cos_theta * (hel.lam1 + hel.lam2) -
              2 * hel.lam3 * kin.cos_theta * (hel.lam1 * hel.lam2 + kin.inv_gamma2));

    if (hel.lam3 == 0)  // longitudinal-transverse
      return 1.5i * G_EM * kin.shat * eft_ext_.s1 * constants::G_F * kin.invA * kin.inv_gamma * kin.sin_theta *
             ((hel.lam2 - hel.lam1) * kin.beta2 - kin.beta * kin.cos_theta * (hel.lam2 + hel.lam1) -
              2 * hel.lam4 * kin.cos_theta * (hel.lam2 * hel.lam1 + kin.inv_gamma2));

    // transverse-transverse
    return 0.75i * G_EM * kin.shat * eft_ext_.s1 * M_SQRT2 * constants::G_F *
           (-kin.inv_gamma2 * kin.beta * (1. + kin.cos_theta2) * (hel.lam1 + hel.lam2) * (hel.lam3 + hel.lam4) +
            2 * kin.sin_theta2 *
                (3. + hel.lam3 * hel.lam4 + hel.lam1 * hel.lam2 * (1 - hel.lam3 * hel.lam4) -
                 kin.beta * (hel.lam1 + hel.lam2) * (hel.lam3 + hel.lam4)) -
            2 * kin.inv_gamma2 *
                (2 + (1 - hel.lam1 * hel.lam2) * hel.lam3 * hel.lam4 -
                 kin.cos_theta2 * (3 + hel.lam1 * hel.lam2 + 2 * hel.lam3 * hel.lam4)));
  }

  std::complex<double> NachtmannAmplitudes::amplitudeWbar(const Kinematics& kin, const Helicities& hel) const {
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return -3 * G_EM * kin.shat * eft_ext_.s1 * M_SQRT2 * constants::G_F * kin.inv_gamma2 * kin.invA *
             kin.sin_theta2 * (hel.lam1 + hel.lam2);

    if (hel.lam4 == 0)  // transverse-longitudinal
      return 1.5 * G_EM * kin.shat * eft_ext_.s1 * constants::G_F * kin.inv_gamma * kin.invA * kin.sin_theta *
             (kin.beta * (hel.lam1 - hel.lam2) * hel.lam3 +
              kin.cos_theta * (2 * kin.beta + (2. - kin.beta2) * (hel.lam1 + hel.lam2) * hel.lam3));

    if (hel.lam3 == 0)  // longitudinal-transverse
      return 1.5 * G_EM * kin.shat * eft_ext_.s1 * constants::G_F * kin.inv_gamma * kin.invA * kin.sin_theta *
             (kin.beta * (hel.lam2 - hel.lam1) * hel.lam4 +
              kin.cos_theta * (2 * kin.beta + (2. - kin.beta2) * (hel.lam2 + hel.lam1) * hel.lam4));

    // transverse-transverse
    return -1.5 * G_EM * kin.shat * eft_ext_.s1 * M_SQRT2 * constants::G_F * kin.invA *
           (2 * kin.sin_theta2 * (hel.lam1 + hel.lam2 - kin.beta * (hel.lam3 + hel.lam4)) +
            kin.inv_gamma2 * ((hel.lam1 + hel.lam2) * (kin.cos_theta2 * (2 + hel.lam3 * hel.lam4) - 1) -
                              kin.beta * (kin.cos_theta2 + hel.lam1 * hel.lam2) * (hel.lam3 + hel.lam4)));
  }

  std::complex<double> NachtmannAmplitudes::amplitudephiW(const Kinematics& kin, const Helicities& hel) const {
    const double invB = 1. / (kin.shat - eft_ext_.mH * eft_ext_.mH);
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return -0.25i * kin.shat * kin.shat * eft_ext_.s1 * eft_ext_.s1 * M_SQRT2 * constants::G_F * invB *
             (1. + kin.beta2) * (1. + hel.lam1 * hel.lam2);

    if (hel.lam4 == 0 || hel.lam3 == 0)  // transverse-longitudinal or longitudinal-transverse
      return 0.;

    // transverse-transverse
    return -0.125i * kin.shat * kin.shat * eft_ext_.s1 * eft_ext_.s1 * M_SQRT2 * constants::G_F * kin.inv_gamma2 *
           invB * (1. + hel.lam1 * hel.lam2) * (1. + hel.lam3 * hel.lam4);
  }

  std::complex<double> NachtmannAmplitudes::amplitudeWB(const Kinematics& kin, const Helicities& hel) const {
    const double invB = 1. / (kin.shat - eft_ext_.mH * eft_ext_.mH);
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return 2.i * G_EM_SQ * kin.invA * eft_ext_.c1() / eft_ext_.s1 *
                 (1 - hel.lam1 * hel.lam2 - 2 * kin.cos_theta2 -
                  kin.gamma2 * (1. + hel.lam1 * hel.lam2) * kin.sin_theta2) +
             0.5 * kin.shat * kin.shat * constants::G_F * M_SQRT2 * invB * eft_ext_.s1 * eft_ext_.c1() *
                 (1. + kin.beta2) * (1. + hel.lam1 * hel.lam2);

    if (hel.lam4 == 0)  // transverse-longitudinal
      return 0.5i * G_EM_SQ * kin.gamma * M_SQRT2 * invB * eft_ext_.c1() / eft_ext_.s1 * kin.sin_theta *
             ((hel.lam2 - hel.lam1) * (1. + kin.inv_gamma2) +
              (kin.beta * float(hel.lam1 + hel.lam2) + 2 * hel.lam3 * (hel.lam1 * hel.lam2 - kin.inv_gamma2)) *
                  kin.cos_theta);

    if (hel.lam3 == 0)  // longitudinal-transverse
      return 0.5i * G_EM_SQ * kin.gamma * M_SQRT2 * invB * eft_ext_.c1() / eft_ext_.s1 * kin.sin_theta *
             ((hel.lam1 - hel.lam2) * (1. + kin.inv_gamma2) +
              (kin.beta * float(hel.lam2 + hel.lam1) + 2 * hel.lam4 * (hel.lam2 * hel.lam1 - kin.inv_gamma2)) *
                  kin.cos_theta);

    // transverse-transverse
    return -0.5i * G_EM_SQ * kin.invA * eft_ext_.c1() / eft_ext_.s1 *
               (kin.beta * float(hel.lam1 + hel.lam2) * (hel.lam3 + hel.lam4) * (1. + kin.cos_theta2) +
                2 * (2 + (hel.lam1 - hel.lam2) * (hel.lam3 - hel.lam4) * kin.cos_theta +
                     ((hel.lam1 * hel.lam2 - 1) * kin.cos_theta2 + 1. + hel.lam1 * hel.lam2) * hel.lam3 * hel.lam4)) +
           0.25i * kin.shat * kin.shat * M_SQRT2 * constants::G_F * kin.inv_gamma2 * invB * eft_ext_.s1 *
               eft_ext_.c1() * (1. + hel.lam1 * hel.lam2) * (1. + hel.lam3 * hel.lam4);
  }

  std::complex<double> NachtmannAmplitudes::amplitudeWbarB(const Kinematics& kin, const Helicities& hel) const {
    CG_WARNING("NachtmannAmplitudes") << "Mode " << mode_ << " is not yet properly handled!";
    const double invB = 1. / (kin.shat - eft_ext_.mH * eft_ext_.mH);
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return 2. * G_EM_SQ * eft_ext_.c1() / eft_ext_.s1 * kin.gamma2 * float(hel.lam1 + hel.lam2) -
             0.5 * kin.shat * kin.shat * M_SQRT2 * constants::G_F /* /e^2 */ * eft_ext_.s1 * eft_ext_.c1() *
                 (1. + kin.beta2) * float(hel.lam1 + hel.lam2);

    if (hel.lam4 == 0)  // transverse-longitudinal
      return 0.5 * G_EM_SQ * kin.invA * kin.gamma * M_SQRT2 * eft_ext_.c1() / eft_ext_.s1 * kin.sin_theta *
             (kin.beta * (hel.lam2 - hel.lam1) * hel.lam3 -
              kin.cos_theta * (2. * kin.beta + kin.beta2 * float(hel.lam1 + hel.lam2) * hel.lam3));

    if (hel.lam3 == 0)  // longitudinal-transverse
      return 0.5 * G_EM_SQ * kin.invA * kin.gamma * M_SQRT2 * eft_ext_.c1() / eft_ext_.s1 * kin.sin_theta *
             (kin.beta * (hel.lam1 - hel.lam2) * hel.lam4 -
              kin.cos_theta * (2. * kin.beta + kin.beta2 * (hel.lam2 + hel.lam1) * hel.lam4));

    // transverse-transverse
    return kin.invA * G_EM_SQ * eft_ext_.c1() * eft_ext_.c1() / eft_ext_.s1 *
               (hel.lam3 * float(hel.lam1 + hel.lam2) + kin.beta * (hel.lam1 * hel.lam2 + kin.cos_theta2)) *
               (hel.lam3 + hel.lam4) -
           0.25 * kin.shat * kin.shat * M_SQRT2 * constants::G_F /* /e^2 */ * kin.inv_gamma2 * invB * eft_ext_.s1 *
               eft_ext_.c1() * eft_ext_.c1() * float(hel.lam1 + hel.lam2) * (1. + hel.lam3 * hel.lam4);
  }

  std::ostream& operator<<(std::ostream& os, const NachtmannAmplitudes::Mode& mode) {
    switch (mode) {
      case NachtmannAmplitudes::Mode::SM:
        return os << "Standard model";
      case NachtmannAmplitudes::Mode::W:
        return os << "W";
      case NachtmannAmplitudes::Mode::Wbar:
        return os << "Wbar";
      case NachtmannAmplitudes::Mode::phiW:
        return os << "phi-W";
      case NachtmannAmplitudes::Mode::phiWbar:
        return os << "phi-Wbar";
      case NachtmannAmplitudes::Mode::phiB:
        return os << "phi-B";
      case NachtmannAmplitudes::Mode::phiBbar:
        return os << "phi-Bbar";
      case NachtmannAmplitudes::Mode::WB:
        return os << "W-B";
      case NachtmannAmplitudes::Mode::WbarB:
        return os << "Wbar-B";
    }
    return os << (int)mode;
  }
}  // namespace cepgen
