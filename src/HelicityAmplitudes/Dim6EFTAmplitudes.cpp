/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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
#include "CepGen/HelicityAmplitudes/FormFactor.h"
#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"
#include "CepGen/Physics/Constants.h"

using namespace cepgen;
using namespace cepgen::amplitudes;
using namespace std::complex_literals;

/// Common contact terms for the dimension-6 anomalous quartic gauge couplings (\f$a_0\f$ and \f$a_C\f$ operators)
class Dim6EFTModel : public HelicityAmplitudes {
public:
  explicit Dim6EFTModel(const ParametersList& params)
      : HelicityAmplitudes(params), form_factor_(FormFactor::build(steer<ParametersList>("formFactor"))) {}

  static ParametersDescription description() {
    auto desc = HelicityAmplitudes::description();
    desc.add("formFactor", FormFactor::description())
        .setDescription("two-parton squared mass form factor of the anomalous couplings");
    return desc;
  }

protected:
  /// Contact amplitude for the (already normalised, in units of the charge squared) \f$a_0\f$ and \f$a_C\f$ couplings
  std::complex<double> eval(double c0, double cc, const Kinematics& kin, const HelicityState& hel) const {
    // Contact terms from the Belanger-Boudjema operators (a_0 and a_C), with
    //   L_0 = -e^2/(16 Lambda^2) c_0 F_{mn}F^{mn} W+^a W-_a
    //   L_C = -e^2/(16 Lambda^2) c_C F_{ma}F^{mb} (W+^a W-_b + W-^a W+_b)
    // i.e. c_0 = 2 a_0^W and c_C = a_C^W for the W bosons (L_0 = -e^2/(8 Lambda^2) a_0^W F_{mn}F^{mn} W+^a W-_a), and
    // c_0 = 2 a_0^Z/cos^2(theta_W) and c_C = a_C^Z/cos^2(theta_W) for the Z bosons (identical particles factor included)
    // amplitude, divided by shat, in a convention-free (real) form
    double val = 0.;
    if (hel[0] == hel[1]) {  // J_z = 0: both a_0 and a_C contribute
      if (hel[2] == 0 && hel[3] == 0)
        val = (1. + kin.beta2) * kin.gamma2 * (c0 / 8. + cc / 16.);
      else if (hel[2] == hel[3])
        val = -(c0 / 8. + cc / 16.);
    } else {  // J_z = +-2: a_C only
      if (hel[2] == 0 && hel[3] == 0)
        val = -kin.gamma2 * kin.sin_theta2 * cc / 16.;
      else if (hel[2] == hel[3])
        val = -kin.sin_theta2 * cc / 32.;
      else if (hel[2] != 0 && hel[3] != 0)
        val = -std::pow(1. + hel[0] * hel[2] * kin.cos_theta, 2) * cc / 32.;
      else if (hel[3] == 0)  // transverse-longitudinal
        val = hel[0] * M_SQRT2 / 32. * kin.gamma * kin.sin_theta * (1. + hel[0] * hel[2] * kin.cos_theta) * cc;
      else  // longitudinal-transverse
        val = hel[0] * M_SQRT2 / 32. * kin.gamma * kin.sin_theta * (1. - hel[0] * hel[3] * kin.cos_theta) * cc;
    }
    // the phase convention (sign flip whenever lam4 == 0) is the one of the SM amplitude above: both are -i times the
    // matrix element of the corresponding Lagrangian terms (checked against explicit polarisation vectors)
    return 1i * (hel[3] == 0 ? -1. : 1.) * kin.shat * val;
  }

protected:
  const std::unique_ptr<FormFactor> form_factor_;
};

/// \f$\gamma\gamma\to W^+W^-\f$ with \f$a_{0,C}^W\f$ quartic couplings (+ SM if requested)
class Wdim6EFTModel final : public Dim6EFTModel {
public:
  explicit Wdim6EFTModel(const ParametersList& params)
      : Dim6EFTModel(params),
        include_sm_(steer<bool>("includeSM")),
        sm_(HelicityAmplitudesFactory::get().build("standardModel")),
        a0w_(steer<double>("a0w")),
        acw_(steer<double>("acw")) {}
  static ParametersDescription description() {
    auto desc = Dim6EFTModel::description();
    desc.setDescription("dimension-6 AQGC (AAWW vertex)");
    desc.add("a0w", 0.).setDescription("a_0^W/Lambda^2 anomalous quartic coupling (in GeV^-2)");
    desc.add("acw", 0.).setDescription("a_C^W/Lambda^2 anomalous quartic coupling (in GeV^-2)");
    desc.add("includeSM", true).setDescription("include the tree-level SM amplitude");
    return desc;
  }
  bool hasAAWW() const override { return true; }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    const auto ff = (*form_factor_)(kin.shat);
    auto ampl = eval(2. * constants::G_EM_SQ * a0w_ * ff, constants::G_EM_SQ * acw_ * ff, kin, hel);
    if (include_sm_)
      ampl += sm_->amplitude(kin, hel);
    return ampl;
  }

private:
  const bool include_sm_;  ///< also include the (tree-level) SM amplitude, for the W pair production
  const std::unique_ptr<HelicityAmplitudes> sm_;  ///< tree-level SM amplitude
  const double a0w_;
  const double acw_;
};
REGISTER_HELICITY_AMPLITUDES("dim6EFTww", Wdim6EFTModel);

/// \f$\gamma\gamma\to ZZ\f$ with \f$a_{0,C}^Z\f$ quartic couplings (no tree-level SM contribution)
class Zdim6EFTModel final : public Dim6EFTModel {
public:
  explicit Zdim6EFTModel(const ParametersList& params)
      : Dim6EFTModel(params),
        a0z_(steer<double>("a0z")),
        acz_(steer<double>("acz")),
        cos2_theta_w_(steer<double>("cos2ThetaW")) {
    if (cos2_theta_w_ <= 0. || cos2_theta_w_ > 1.)
      throw CG_FATAL("Zdim6EFTModel") << "Invalid squared cosine of the weak mixing angle: " << cos2_theta_w_ << ".";
  }
  static ParametersDescription description() {
    auto desc = Dim6EFTModel::description();
    desc.setDescription("dimension-6 AQGC (AAZZ vertex)");
    desc.add("a0z", 0.).setDescription("a_0^Z/Lambda^2 anomalous quartic coupling (in GeV^-2)");
    desc.add("acz", 0.).setDescription("a_C^Z/Lambda^2 anomalous quartic coupling (in GeV^-2)");
    desc.add("cos2ThetaW", 0.76878).setDescription("cos^2(theta_W), entering the Z quartic couplings normalisation");
    return desc;
  }
  bool hasAAZZ() const override { return true; }
  std::complex<double> amplitude(const Kinematics& kin, const HelicityState& hel) const override {
    const auto ff = (*form_factor_)(kin.shat);
    return eval(constants::G_EM_SQ * a0z_ * ff / cos2_theta_w_ * 2.,  // factor 2 from identical Z bosons
                constants::G_EM_SQ * acz_ * ff / cos2_theta_w_,       // a_C already symmetrised
                kin,
                hel);
  }

private:
  const double a0z_;
  const double acz_;
  const double cos2_theta_w_;  ///< squared cosine of the weak mixing angle, used in the Z couplings
};
REGISTER_HELICITY_AMPLITUDES("dim6EFTzz", Zdim6EFTModel);
