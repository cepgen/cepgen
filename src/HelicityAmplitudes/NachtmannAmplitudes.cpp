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

#include <algorithm>

#include "CepGen/Core/Exception.h"
#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"
#include "CepGen/Physics/Constants.h"

using namespace cepgen;
using namespace cepgen::amplitudes;
using namespace std::complex_literals;

///////////////////////////////////////////////////////////////////////////////////////////////////

struct StandardModel : HelicityAmplitudes {
  explicit StandardModel(const ParametersList& params) : HelicityAmplitudes(params) {}
  static ParametersDescription description() {
    auto desc = HelicityAmplitudes::description();
    desc.setDescription("Standard model");
    return desc;
  }
  bool hasAAWW() const override { return true; }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    if (hel[2] == 0 && hel[3] == 0)  // longitudinal-longitudinal
      return 1i * constants::G_EM_SQ * kin.invA * kin.inv_gamma2 *
             ((kin.gamma2 + 1.) * (1. - hel[0] * hel[1]) * kin.sin_theta2 - (1. + hel[0] * hel[1]));

    if (hel[3] == 0)  // transverse-longitudinal
      return -1i * constants::G_EM_SQ * M_SQRT2 * kin.invA * kin.inv_gamma * static_cast<double>(hel[0] - hel[1]) *
             (1. + hel[0] * hel[2] * kin.cos_theta) * kin.sin_theta;

    if (hel[2] == 0)  // longitudinal-transverse
      return amplitude(kin, {hel[1], hel[0], hel[3], hel[2]});

    // transverse-transverse
    return -0.5i * constants::G_EM_SQ * kin.invA *
           (2. * kin.beta * static_cast<double>(hel[0] + hel[1]) * (hel[2] + hel[3]) -
            kin.inv_gamma2 * (1. + hel[2] * hel[3]) * (2. * hel[0] * hel[1] + (1. - hel[0] * hel[1]) * kin.cos_theta2) +
            (1. + hel[0] * hel[1] * hel[2] * hel[3]) * (3. + hel[0] * hel[1]) +
            2. * (hel[0] - hel[1]) * (hel[2] - hel[3]) * kin.cos_theta +
            (1. - hel[0] * hel[1]) * (1. - hel[2] * hel[3]) * kin.cos_theta2);
  }
};
REGISTER_HELICITY_AMPLITUDES("standardModel", StandardModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

class WModel : public HelicityAmplitudes {
public:
  explicit WModel(const ParametersList& params) : HelicityAmplitudes(params) {}
  static ParametersDescription description() {
    auto desc = HelicityAmplitudes::description();
    desc.setDescription("W");
    return desc;
  }
  bool hasAAWW() const override { return true; }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    if (hel[2] == 0 && hel[3] == 0)  // longitudinal-longitudinal
      return 3i * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * M_SQRT2 * constants::G_F * kin.invA *
             kin.inv_gamma2 * kin.sin_theta2 * (1. + hel[0] * hel[1]);

    if (hel[3] == 0)  // transverse-longitudinal
      return 1.5i * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * constants::G_F * kin.invA * kin.inv_gamma *
             kin.sin_theta *
             ((hel[0] - hel[1]) * kin.beta2 - kin.beta * kin.cos_theta * (hel[0] + hel[1]) -
              2 * hel[2] * kin.cos_theta * (hel[0] * hel[1] + kin.inv_gamma2));

    if (hel[2] == 0)  // longitudinal-transverse
      return amplitude(kin, {hel[1], hel[0], hel[3], hel[2]});

    // transverse-transverse
    return 0.75i * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * M_SQRT2 * constants::G_F *
           (-kin.inv_gamma2 * kin.beta * (1. + kin.cos_theta2) * (hel[0] + hel[1]) * (hel[2] + hel[3]) +
            2 * kin.sin_theta2 *
                (3. + hel[2] * hel[3] + hel[0] * hel[1] * (1 - hel[2] * hel[3]) -
                 kin.beta * (hel[0] + hel[1]) * (hel[2] + hel[3])) -
            2 * kin.inv_gamma2 *
                (2 + (1 - hel[0] * hel[1]) * hel[2] * hel[3] -
                 kin.cos_theta2 * (3 + hel[0] * hel[1] + 2 * hel[2] * hel[3])));
  }
};
REGISTER_HELICITY_AMPLITUDES("W", WModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

class WbarModel : public HelicityAmplitudes {
public:
  explicit WbarModel(const ParametersList& params) : HelicityAmplitudes(params) {}
  static ParametersDescription description() {
    auto desc = HelicityAmplitudes::description();
    desc.setDescription("Wbar");
    return desc;
  }
  bool hasAAWW() const override { return true; }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    if (hel[2] == 0 && hel[3] == 0)  // longitudinal-longitudinal
      return {-3 * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * M_SQRT2 * constants::G_F * kin.inv_gamma2 *
                  kin.invA * kin.sin_theta2 * (hel[0] + hel[1]),
              0.};

    if (hel[3] == 0)  // transverse-longitudinal
      return {1.5 * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * constants::G_F * kin.inv_gamma * kin.invA *
                  kin.sin_theta *
                  (kin.beta * (hel[0] - hel[1]) * hel[2] +
                   kin.cos_theta * (2 * kin.beta + (2. - kin.beta2) * (hel[0] + hel[1]) * hel[2])),
              0.};

    if (hel[2] == 0)  // longitudinal-transverse
      return amplitude(kin, {hel[1], hel[0], hel[3], hel[2]});

    // transverse-transverse
    return {-1.5 * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * M_SQRT2 * constants::G_F * kin.invA *
                (2 * kin.sin_theta2 * (hel[0] + hel[1] - kin.beta * (hel[2] + hel[3])) +
                 kin.inv_gamma2 * ((hel[0] + hel[1]) * (kin.cos_theta2 * (2 + hel[2] * hel[3]) - 1) -
                                   kin.beta * (kin.cos_theta2 + hel[0] * hel[1]) * (hel[2] + hel[3]))),
            0.};
  }
};
REGISTER_HELICITY_AMPLITUDES("Wbar", WbarModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

class PhiWModel : public HelicityAmplitudes {
public:
  explicit PhiWModel(const ParametersList& params) : HelicityAmplitudes(params), m_h_(steer<double>("mH")) {}
  static ParametersDescription description() {
    auto desc = HelicityAmplitudes::description();
    desc.setDescription("phi-W");
    desc.add("mH", 0.).setDescription("Higgs mass (in GeV/c2)");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    const double invB = 1. / (kin.shat - m_h_ * m_h_);
    if (hel[2] == 0 && hel[3] == 0)  // longitudinal-longitudinal
      return -0.25i * kin.shat2 * s1() * s1() * M_SQRT2 * constants::G_F * invB * (1. + kin.beta2) *
             (1. + hel[0] * hel[1]);

    if (hel[3] == 0 || hel[2] == 0)  // transverse-longitudinal or longitudinal-transverse
      return {0., 0.};

    // transverse-transverse
    return -0.125i * kin.shat2 * s1() * s1() * M_SQRT2 * constants::G_F * kin.inv_gamma2 * invB *
           (1. + hel[0] * hel[1]) * (1. + hel[2] * hel[3]);
  }

private:
  const double m_h_;
};
REGISTER_HELICITY_AMPLITUDES("PhiW", PhiWModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

struct PhiWbarModel : PhiWModel {
  using PhiWModel::PhiWModel;
  static ParametersDescription description() {
    auto desc = PhiWModel::description();
    desc.setDescription("phi-Wbar");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    return 2i * static_cast<double>(hel[0]) * PhiWModel::amplitude(kin, hel);
  }
};
REGISTER_HELICITY_AMPLITUDES("PhiWbar", PhiWbarModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

struct PhiBModel : PhiWModel {
  explicit PhiBModel(const ParametersList& params) : PhiWModel(params) {}
  static ParametersDescription description() {
    auto desc = PhiWModel::description();
    desc.setDescription("phi-B");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    return std::pow(c1() / s1(), 2) * PhiWModel::amplitude(kin, hel);
  }
};
REGISTER_HELICITY_AMPLITUDES("PhiB", PhiBModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

struct PhiBbarModel : PhiWModel {
  explicit PhiBbarModel(const ParametersList& params) : PhiWModel(params) {}
  static ParametersDescription description() {
    auto desc = PhiWModel::description();
    desc.setDescription("phi-Bbar");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    return 2i * static_cast<double>(hel[0]) * std::pow(c1() / s1(), 2) * PhiWModel::amplitude(kin, hel);
  }
};
REGISTER_HELICITY_AMPLITUDES("PhiBbar", PhiBbarModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

class WBModel : public HelicityAmplitudes {
public:
  explicit WBModel(const ParametersList& params) : HelicityAmplitudes(params), m_h_(steer<double>("mH")) {}
  static ParametersDescription description() {
    auto desc = HelicityAmplitudes::description();
    desc.setDescription("W-b");
    desc.add("mH", 0.).setDescription("Higgs mass (in GeV/c2)");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    const auto invB = 1. / (kin.shat - m_h_ * m_h_);
    if (hel[2] == 0 && hel[3] == 0)  // longitudinal-longitudinal
      return 2i * constants::G_EM_SQ * kin.invA * c1() / s1() *
                 (1 - hel[0] * hel[1] - 2 * kin.cos_theta2 - kin.gamma2 * (1. + hel[0] * hel[1]) * kin.sin_theta2) +
             0.5i * kin.shat2 * constants::G_F * M_SQRT2 * invB * s1() * c1() * (1. + kin.beta2) *
                 (1. + hel[0] * hel[1]);

    if (hel[3] == 0)  // transverse-longitudinal
      return 0.5i * constants::G_EM_SQ * kin.gamma * M_SQRT2 * kin.invA * c1() / s1() * kin.sin_theta *
             ((hel[1] - hel[0]) * (1. + kin.inv_gamma2) +
              (kin.beta * static_cast<double>(hel[0] + hel[1]) + 2 * hel[2] * (hel[0] * hel[1] - kin.inv_gamma2)) *
                  kin.cos_theta);

    if (hel[2] == 0)  // longitudinal-transverse
      return amplitude(kin, {hel[1], hel[0], hel[3], hel[2]});

    // transverse-transverse
    return -0.5i * constants::G_EM_SQ * kin.invA * c1() / s1() *
               (kin.beta * static_cast<double>(hel[0] + hel[1]) * (hel[2] + hel[3]) * (1. + kin.cos_theta2) +
                2 * (2 + (hel[0] - hel[1]) * (hel[2] - hel[3]) * kin.cos_theta +
                     ((hel[0] * hel[1] - 1) * kin.cos_theta2 + 1. + hel[0] * hel[1]) * hel[2] * hel[3])) +
           0.25i * kin.shat2 * M_SQRT2 * constants::G_F * kin.inv_gamma2 * invB * s1() * c1() * (1. + hel[0] * hel[1]) *
               (1. + hel[2] * hel[3]);
  }

private:
  const double m_h_;
};
REGISTER_HELICITY_AMPLITUDES("WB", WBModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

class WbarBModel : public HelicityAmplitudes {
public:
  explicit WbarBModel(const ParametersList& params) : HelicityAmplitudes(params), m_h_(steer<double>("mH")) {
    CG_WARNING("HelicityAmplitudes") << "Mode Wbar-B is not yet properly handled!";
  }
  static ParametersDescription description() {
    auto desc = HelicityAmplitudes::description();
    desc.setDescription("Wbar-B");
    desc.add("mH", 0.).setDescription("Higgs mass (in GeV/c2)");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    const auto invB = 1. / (kin.shat - m_h_ * m_h_);
    if (hel[2] == 0 && hel[3] == 0)  // longitudinal-longitudinal
      return {2. * constants::G_EM_SQ * c1() / s1() * kin.gamma2 * static_cast<double>(hel[0] + hel[1]) -
                  0.5 * kin.shat2 * M_SQRT2 * constants::G_F /* /e^2 */ * invB * s1() * c1() * (1. + kin.beta2) *
                      static_cast<double>(hel[0] + hel[1]),
              0.};

    if (hel[3] == 0)  // transverse-longitudinal
      return {0.5 * constants::G_EM_SQ * kin.invA * kin.gamma * M_SQRT2 * c1() / s1() * kin.sin_theta *
                  (kin.beta * (hel[1] - hel[0]) * hel[2] -
                   kin.cos_theta * (2. * kin.beta + kin.beta2 * static_cast<double>(hel[0] + hel[1]) * hel[2])),
              0.};

    if (hel[2] == 0)  // longitudinal-transverse
      return amplitude(kin, {hel[1], hel[0], hel[3], hel[2]});

    // transverse-transverse
    return {kin.invA * constants::G_EM_SQ * c1() * c1() / s1() *
                    (hel[2] * static_cast<double>(hel[0] + hel[1]) + kin.beta * (hel[0] * hel[1] + kin.cos_theta2)) *
                    (hel[2] + hel[3]) -
                0.25 * kin.shat2 * M_SQRT2 * constants::G_F /* /e^2 */ * kin.inv_gamma2 * invB * s1() * c1() * c1() *
                    static_cast<double>(hel[0] + hel[1]) * (1. + hel[2] * hel[3]),
            0.};
  }

private:
  const double m_h_;
};
REGISTER_HELICITY_AMPLITUDES("WbarB", WbarBModel);
