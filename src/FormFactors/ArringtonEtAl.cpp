/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2025  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"

namespace cepgen::formfac {
  /// \cite Arrington:2007ux
  class ArringtonEtAl final : public Parameterisation {
  public:
    explicit ArringtonEtAl(const ParametersList& params) : Parameterisation(params), mode_(steer<int>("mode")) {
      switch (mode_) {
        case 0:  // original
          a_e_ = {3.439, -1.602, 0.068};
          b_e_ = {15.055, 48.061, 99.304, 0.012, 8.650};
          a_m_ = {-1.465, 1.260, 0.262};
          b_m_ = {9.627, 0., 0., 11.179, 13.245};
          break;
        case 1:  // fit of quoted Ge+dGe values
          a_e_ = {4.309, -1.108, -0.324};
          b_e_ = {15.340, 58.321, 124.11, 3.927, 0.589};
          a_m_ = {-1.472, 1.210, 0.334};
          b_m_ = {9.486, 0., 0., 9.440, 15.416};
          break;
        case 2:  // fit of quoted Ge-dGe values
          a_e_ = {4.286, -1.281, -0.486};
          b_e_ = {16.308, 54.535, 138.03, 7.005, 0.014};
          a_m_ = {-1.374, 1.080, 0.124};
          b_m_ = {10.003, 0., 0., 7.680, 9.009};
          break;
        case 3:  // fit of quoted Ge values
          a_e_ = {4.109, -1.052, -0.375};
          b_e_ = {15.602, 55.519, 123.96, 11.403, 1.931};
          a_m_ = {-1.436, 1.196, 0.210};
          b_m_ = {9.721, 0., 0., 9.623, 11.817};
          break;
        default:
          throw CG_FATAL("ArringtonEtAl") << "Invalid parameterisation mode: " << mode_ << ".";
      }
    }

    inline static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Arrington et al.");
      desc.add<int>("mode", 0)
          .setDescription("Parameterisation mode")
          .allow(0, "original")
          .allow(1, "fit of quoted Ge+dGe values")
          .allow(2, "fit of quoted Ge-dGe values")
          .allow(3, "fit of quoted Ge values");
      return desc;
    }

  private:
    inline void eval() override {
      const double tau_val = tau(q2_);

      double num_e = 1., den_e = 1.;
      for (size_t i = 0; i < a_e_.size(); ++i)
        num_e += a_e_.at(i) * std::pow(tau_val, 1. + i);
      for (size_t i = 0; i < b_e_.size(); ++i)
        den_e += b_e_.at(i) * std::pow(tau_val, 1. + i);
      const auto ge = num_e / den_e;

      double num_m = 1., den_m = 1.;
      for (size_t i = 0; i < a_m_.size(); ++i)
        num_m += a_m_.at(i) * std::pow(tau_val, 1. + i);
      for (size_t i = 0; i < b_m_.size(); ++i)
        den_m += b_m_.at(i) * std::pow(tau_val, 1. + i);
      const auto gm = MU * num_m / den_m;

      setGEGM(ge, gm);
    }

    const int mode_;
    std::vector<double> a_e_, b_e_;
    std::vector<double> a_m_, b_m_;
  };
}  // namespace cepgen::formfac
using cepgen::formfac::ArringtonEtAl;
REGISTER_FORMFACTORS("Arrington", ArringtonEtAl);
