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
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/NachtmannAmplitudes.h"

using namespace cepgen;
using namespace cepgen::amplitudes;
using namespace std::complex_literals;

///////////////////////////////////////////////////////////////////////////////////////////////////

Kinematics::Kinematics(double mv2, double s_hat, double t_hat, double u_hat)
    : shat(s_hat),
      that(t_hat),
      uhat(u_hat),
      mv2_(mv2),
      shat2(s_hat * s_hat),
      beta2(1. - 4. * mv2_ / s_hat),
      beta(sqrt(beta2)),
      inv_gamma2(1. - beta2),
      gamma2(1. / inv_gamma2),
      gamma(sqrt(gamma2)),
      inv_gamma(1. / gamma) {
  setCosTheta((t_hat - u_hat) / s_hat / beta);
}

Kinematics Kinematics::fromSCosTheta(double shat, double cos_theta, double mv2) {
  Kinematics kin(mv2, shat, 0., 0.);
  kin.setCosTheta(cos_theta);
  return kin;
}

void Kinematics::setCosTheta(double cth) {
  // with off-shell photons, (t-u)/(s*beta) may exceed unity, possibly reaching the 1/A pole at beta*|cos(theta)|=1
  cos_theta = std::clamp(cth, -1., 1.);
  cos_theta2 = cos_theta * cos_theta;
  sin_theta2 = 1. - cos_theta2;
  sin_theta = sqrt(fabs(sin_theta2));
  invA = 1. / (1. - beta2 * cos_theta2);
}

bool Kinematics::operator!=(const Kinematics& oth) const {
  // checking only the base variables as all others are computed from these three
  return shat != oth.shat || that != oth.that || uhat != oth.uhat;
}

namespace cepgen::amplitudes {
  std::ostream& operator<<(std::ostream& os, const Kinematics& kin) {
    return os << "Kin{mV2=" << kin.mv2_ << ",shat=" << kin.shat << ",that=" << kin.that << ",uhat=" << kin.uhat
              << ",beta=" << kin.beta << ",gamma=" << kin.gamma << ",cos(theta)=" << kin.cos_theta
              << "->1/A=" << kin.invA << "}";
  }
}  // namespace cepgen::amplitudes

///////////////////////////////////////////////////////////////////////////////////////////////////

NachtmannAmplitudes::NachtmannAmplitudes(const ParametersList& params)
    : NamedModule(params), s1_(steer<double>("s1")) {}

std::complex<double> NachtmannAmplitudes::operator()(
    const Kinematics& kin, short lam1, short lam2, short lam3, short lam4) const {
  return amplitude(kin, HelicityStates{lam1, lam2, lam3, lam4});
}

ParametersDescription NachtmannAmplitudes::description() {
  auto desc = ParametersDescription();
  desc.setDescription("generic helicity model");
  desc.add("s1", 0.);
  return desc;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

struct StandardModel : NachtmannAmplitudes {
  explicit StandardModel(const ParametersList& params) : NachtmannAmplitudes(params) {}
  static ParametersDescription description() {
    auto desc = NachtmannAmplitudes::description();
    desc.setDescription("Standard model");
    return desc;
  }
  bool hasAAWW() const override { return true; }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityStates& hel) const override {
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return 1i * constants::G_EM_SQ * kin.invA * kin.inv_gamma2 *
             ((kin.gamma2 + 1.) * (1. - hel.lam1 * hel.lam2) * kin.sin_theta2 - (1. + hel.lam1 * hel.lam2));

    if (hel.lam4 == 0)  // transverse-longitudinal
      return -1i * constants::G_EM_SQ * M_SQRT2 * kin.invA * kin.inv_gamma * static_cast<double>(hel.lam1 - hel.lam2) *
             (1. + hel.lam1 * hel.lam3 * kin.cos_theta) * kin.sin_theta;

    if (hel.lam3 == 0)  // longitudinal-transverse
      return amplitude(kin, {hel.lam2, hel.lam1, hel.lam4, hel.lam3});

    // transverse-transverse
    return -0.5i * constants::G_EM_SQ * kin.invA *
           (2. * kin.beta * static_cast<double>(hel.lam1 + hel.lam2) * (hel.lam3 + hel.lam4) -
            kin.inv_gamma2 * (1. + hel.lam3 * hel.lam4) *
                (2. * hel.lam1 * hel.lam2 + (1. - hel.lam1 * hel.lam2) * kin.cos_theta2) +
            (1. + hel.lam1 * hel.lam2 * hel.lam3 * hel.lam4) * (3. + hel.lam1 * hel.lam2) +
            2. * (hel.lam1 - hel.lam2) * (hel.lam3 - hel.lam4) * kin.cos_theta +
            (1. - hel.lam1 * hel.lam2) * (1. - hel.lam3 * hel.lam4) * kin.cos_theta2);
  }
};
REGISTER_NACHTMANN_AMPLITUDES("standardModel", StandardModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

class WModel : public NachtmannAmplitudes {
public:
  explicit WModel(const ParametersList& params) : NachtmannAmplitudes(params) {}
  static ParametersDescription description() {
    auto desc = NachtmannAmplitudes::description();
    desc.setDescription("W");
    return desc;
  }
  bool hasAAWW() const override { return true; }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityStates& hel) const override {
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return 3i * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * M_SQRT2 * constants::G_F * kin.invA *
             kin.inv_gamma2 * kin.sin_theta2 * (1. + hel.lam1 * hel.lam2);

    if (hel.lam4 == 0)  // transverse-longitudinal
      return 1.5i * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * constants::G_F * kin.invA * kin.inv_gamma *
             kin.sin_theta *
             ((hel.lam1 - hel.lam2) * kin.beta2 - kin.beta * kin.cos_theta * (hel.lam1 + hel.lam2) -
              2 * hel.lam3 * kin.cos_theta * (hel.lam1 * hel.lam2 + kin.inv_gamma2));

    if (hel.lam3 == 0)  // longitudinal-transverse
      return amplitude(kin, {hel.lam2, hel.lam1, hel.lam4, hel.lam3});

    // transverse-transverse
    return 0.75i * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * M_SQRT2 * constants::G_F *
           (-kin.inv_gamma2 * kin.beta * (1. + kin.cos_theta2) * (hel.lam1 + hel.lam2) * (hel.lam3 + hel.lam4) +
            2 * kin.sin_theta2 *
                (3. + hel.lam3 * hel.lam4 + hel.lam1 * hel.lam2 * (1 - hel.lam3 * hel.lam4) -
                 kin.beta * (hel.lam1 + hel.lam2) * (hel.lam3 + hel.lam4)) -
            2 * kin.inv_gamma2 *
                (2 + (1 - hel.lam1 * hel.lam2) * hel.lam3 * hel.lam4 -
                 kin.cos_theta2 * (3 + hel.lam1 * hel.lam2 + 2 * hel.lam3 * hel.lam4)));
  }
};
REGISTER_NACHTMANN_AMPLITUDES("W", WModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

class WbarModel : public NachtmannAmplitudes {
public:
  explicit WbarModel(const ParametersList& params) : NachtmannAmplitudes(params) {}
  static ParametersDescription description() {
    auto desc = NachtmannAmplitudes::description();
    desc.setDescription("Wbar");
    return desc;
  }
  bool hasAAWW() const override { return true; }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityStates& hel) const override {
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return {-3 * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * M_SQRT2 * constants::G_F * kin.inv_gamma2 *
                  kin.invA * kin.sin_theta2 * (hel.lam1 + hel.lam2),
              0.};

    if (hel.lam4 == 0)  // transverse-longitudinal
      return {1.5 * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * constants::G_F * kin.inv_gamma * kin.invA *
                  kin.sin_theta *
                  (kin.beta * (hel.lam1 - hel.lam2) * hel.lam3 +
                   kin.cos_theta * (2 * kin.beta + (2. - kin.beta2) * (hel.lam1 + hel.lam2) * hel.lam3)),
              0.};

    if (hel.lam3 == 0)  // longitudinal-transverse
      return amplitude(kin, {hel.lam2, hel.lam1, hel.lam4, hel.lam3});

    // transverse-transverse
    return {-1.5 * std::sqrt(constants::G_EM_SQ) * kin.shat * s1() * M_SQRT2 * constants::G_F * kin.invA *
                (2 * kin.sin_theta2 * (hel.lam1 + hel.lam2 - kin.beta * (hel.lam3 + hel.lam4)) +
                 kin.inv_gamma2 * ((hel.lam1 + hel.lam2) * (kin.cos_theta2 * (2 + hel.lam3 * hel.lam4) - 1) -
                                   kin.beta * (kin.cos_theta2 + hel.lam1 * hel.lam2) * (hel.lam3 + hel.lam4))),
            0.};
  }
};
REGISTER_NACHTMANN_AMPLITUDES("Wbar", WbarModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

class PhiWModel : public NachtmannAmplitudes {
public:
  explicit PhiWModel(const ParametersList& params) : NachtmannAmplitudes(params), m_h_(steer<double>("mH")) {}
  static ParametersDescription description() {
    auto desc = NachtmannAmplitudes::description();
    desc.setDescription("phi-W");
    desc.add("mH", 0.).setDescription("Higgs mass (in GeV/c2)");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityStates& hel) const override {
    const double invB = 1. / (kin.shat - m_h_ * m_h_);
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return -0.25i * kin.shat2 * s1() * s1() * M_SQRT2 * constants::G_F * invB * (1. + kin.beta2) *
             (1. + hel.lam1 * hel.lam2);

    if (hel.lam4 == 0 || hel.lam3 == 0)  // transverse-longitudinal or longitudinal-transverse
      return {0., 0.};

    // transverse-transverse
    return -0.125i * kin.shat2 * s1() * s1() * M_SQRT2 * constants::G_F * kin.inv_gamma2 * invB *
           (1. + hel.lam1 * hel.lam2) * (1. + hel.lam3 * hel.lam4);
  }

private:
  const double m_h_;
};
REGISTER_NACHTMANN_AMPLITUDES("PhiW", PhiWModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

struct PhiWbarModel : PhiWModel {
  using PhiWModel::PhiWModel;
  static ParametersDescription description() {
    auto desc = PhiWModel::description();
    desc.setDescription("phi-Wbar");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityStates& hel) const override {
    return 2i * static_cast<double>(hel.lam1) * PhiWModel::amplitude(kin, hel);
  }
};
REGISTER_NACHTMANN_AMPLITUDES("PhiWbar", PhiWbarModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

struct PhiBModel : PhiWModel {
  explicit PhiBModel(const ParametersList& params) : PhiWModel(params) {}
  static ParametersDescription description() {
    auto desc = PhiWModel::description();
    desc.setDescription("phi-B");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityStates& hel) const override {
    return std::pow(c1() / s1(), 2) * PhiWModel::amplitude(kin, hel);
  }
};
REGISTER_NACHTMANN_AMPLITUDES("PhiB", PhiBModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

struct PhiBbarModel : PhiWModel {
  explicit PhiBbarModel(const ParametersList& params) : PhiWModel(params) {}
  static ParametersDescription description() {
    auto desc = PhiWModel::description();
    desc.setDescription("phi-Bbar");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityStates& hel) const override {
    return 2i * static_cast<double>(hel.lam1) * std::pow(c1() / s1(), 2) * PhiWModel::amplitude(kin, hel);
  }
};
REGISTER_NACHTMANN_AMPLITUDES("PhiBbar", PhiBbarModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

class WBModel : public NachtmannAmplitudes {
public:
  explicit WBModel(const ParametersList& params) : NachtmannAmplitudes(params), m_h_(steer<double>("mH")) {}
  static ParametersDescription description() {
    auto desc = NachtmannAmplitudes::description();
    desc.setDescription("W-b");
    desc.add("mH", 0.).setDescription("Higgs mass (in GeV/c2)");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityStates& hel) const override {
    const auto invB = 1. / (kin.shat - m_h_ * m_h_);
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return 2i * constants::G_EM_SQ * kin.invA * c1() / s1() *
                 (1 - hel.lam1 * hel.lam2 - 2 * kin.cos_theta2 -
                  kin.gamma2 * (1. + hel.lam1 * hel.lam2) * kin.sin_theta2) +
             0.5i * kin.shat2 * constants::G_F * M_SQRT2 * invB * s1() * c1() * (1. + kin.beta2) *
                 (1. + hel.lam1 * hel.lam2);

    if (hel.lam4 == 0)  // transverse-longitudinal
      return 0.5i * constants::G_EM_SQ * kin.gamma * M_SQRT2 * kin.invA * c1() / s1() * kin.sin_theta *
             ((hel.lam2 - hel.lam1) * (1. + kin.inv_gamma2) + (kin.beta * static_cast<double>(hel.lam1 + hel.lam2) +
                                                               2 * hel.lam3 * (hel.lam1 * hel.lam2 - kin.inv_gamma2)) *
                                                                  kin.cos_theta);

    if (hel.lam3 == 0)  // longitudinal-transverse
      return amplitude(kin, {hel.lam2, hel.lam1, hel.lam4, hel.lam3});

    // transverse-transverse
    return -0.5i * constants::G_EM_SQ * kin.invA * c1() / s1() *
               (kin.beta * static_cast<double>(hel.lam1 + hel.lam2) * (hel.lam3 + hel.lam4) * (1. + kin.cos_theta2) +
                2 * (2 + (hel.lam1 - hel.lam2) * (hel.lam3 - hel.lam4) * kin.cos_theta +
                     ((hel.lam1 * hel.lam2 - 1) * kin.cos_theta2 + 1. + hel.lam1 * hel.lam2) * hel.lam3 * hel.lam4)) +
           0.25i * kin.shat2 * M_SQRT2 * constants::G_F * kin.inv_gamma2 * invB * s1() * c1() *
               (1. + hel.lam1 * hel.lam2) * (1. + hel.lam3 * hel.lam4);
  }

private:
  const double m_h_;
};
REGISTER_NACHTMANN_AMPLITUDES("WB", WBModel);

///////////////////////////////////////////////////////////////////////////////////////////////////

class WbarBModel : public NachtmannAmplitudes {
public:
  explicit WbarBModel(const ParametersList& params) : NachtmannAmplitudes(params), m_h_(steer<double>("mH")) {
    CG_WARNING("NachtmannAmplitudes") << "Mode Wbar-B is not yet properly handled!";
  }
  static ParametersDescription description() {
    auto desc = NachtmannAmplitudes::description();
    desc.setDescription("Wbar-B");
    desc.add("mH", 0.).setDescription("Higgs mass (in GeV/c2)");
    return desc;
  }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityStates& hel) const override {
    const auto invB = 1. / (kin.shat - m_h_ * m_h_);
    if (hel.lam3 == 0 && hel.lam4 == 0)  // longitudinal-longitudinal
      return {2. * constants::G_EM_SQ * c1() / s1() * kin.gamma2 * static_cast<double>(hel.lam1 + hel.lam2) -
                  0.5 * kin.shat2 * M_SQRT2 * constants::G_F /* /e^2 */ * invB * s1() * c1() * (1. + kin.beta2) *
                      static_cast<double>(hel.lam1 + hel.lam2),
              0.};

    if (hel.lam4 == 0)  // transverse-longitudinal
      return {0.5 * constants::G_EM_SQ * kin.invA * kin.gamma * M_SQRT2 * c1() / s1() * kin.sin_theta *
                  (kin.beta * (hel.lam2 - hel.lam1) * hel.lam3 -
                   kin.cos_theta * (2. * kin.beta + kin.beta2 * static_cast<double>(hel.lam1 + hel.lam2) * hel.lam3)),
              0.};

    if (hel.lam3 == 0)  // longitudinal-transverse
      return amplitude(kin, {hel.lam2, hel.lam1, hel.lam4, hel.lam3});

    // transverse-transverse
    return {kin.invA * constants::G_EM_SQ * c1() * c1() / s1() *
                    (hel.lam3 * static_cast<double>(hel.lam1 + hel.lam2) +
                     kin.beta * (hel.lam1 * hel.lam2 + kin.cos_theta2)) *
                    (hel.lam3 + hel.lam4) -
                0.25 * kin.shat2 * M_SQRT2 * constants::G_F /* /e^2 */ * kin.inv_gamma2 * invB * s1() * c1() * c1() *
                    static_cast<double>(hel.lam1 + hel.lam2) * (1. + hel.lam3 * hel.lam4),
            0.};
  }

private:
  const double m_h_;
};
REGISTER_NACHTMANN_AMPLITUDES("WbarB", WbarBModel);
