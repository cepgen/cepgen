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
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"

namespace cepgen::formfac {
  /// \cite Mergell:1995bf
  class MergellEtAl final : public Parameterisation {
  public:
    explicit MergellEtAl(const ParametersList& params)
        : Parameterisation(params),
          rho1_coefficients_(steer<std::vector<double> >("rho1Coefficients")),
          rho2_coefficients_(steer<std::vector<double> >("rho2Coefficients")),
          d1_terms_(steer<std::vector<double> >("d1terms")),
          d2_terms_(steer<std::vector<double> >("d2terms")),
          fs1_coefficients_(steer<std::vector<double> >("fs1Coefficients")),
          fs2_coefficients_(steer<std::vector<double> >("fs2Coefficients")),
          fv1_coefficients_(steer<std::vector<double> >("fv1Coefficients")),
          fv2_coefficients_(steer<std::vector<double> >("fv2Coefficients")),
          inv_q20_(steer<double>("q20inv")),
          lambda_sq_(steer<double>("LambdaSq")),
          gamma_(steer<double>("gamma")) {
      if (rho1_coefficients_.size() < 4)
        throw CG_FATAL("MergellEtAl") << "Invalid multiplicity of rho1 coefficients. Got " << rho1_coefficients_.size()
                                      << ", should have 4.";
      if (rho2_coefficients_.size() < 4)
        throw CG_FATAL("MergellEtAl") << "Invalid multiplicity of rho1 coefficients. Got " << rho1_coefficients_.size()
                                      << ", should have 4.";
      if (d1_terms_.size() < 3)
        throw CG_FATAL("MergellEtAl") << "Invalid multiplicity of d1 terms. Got " << d1_terms_.size()
                                      << ", should have 3.";
      if (d2_terms_.size() < 3)
        throw CG_FATAL("MergellEtAl") << "Invalid multiplicity of d2 terms. Got " << d1_terms_.size()
                                      << ", should have 3.";
      if (fs1_coefficients_.size() < 3)
        throw CG_FATAL("MergellEtAl") << "Invalid multiplicity of Fs1 coefficients. Got " << fs1_coefficients_.size()
                                      << ", should have 3.";
      if (fs2_coefficients_.size() < 3)
        throw CG_FATAL("MergellEtAl") << "Invalid multiplicity of Fs2 coefficients. Got " << fs2_coefficients_.size()
                                      << ", should have 3.";
      if (fv1_coefficients_.size() < 3)
        throw CG_FATAL("MergellEtAl") << "Invalid multiplicity of Fv1 coefficients. Got " << fv1_coefficients_.size()
                                      << ", should have 3.";
      if (fv2_coefficients_.size() < 3)
        throw CG_FATAL("MergellEtAl") << "Invalid multiplicity of Fv2 coefficients. Got " << fv2_coefficients_.size()
                                      << ", should have 3.";
    }

    inline static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Mergell et al.");
      desc.add("q20inv", 1. / 0.35);
      desc.add("rho1Coefficients", std::vector{1.0317, 0.0875, 0.3176, 0.5496});
      desc.add("rho2Coefficients", std::vector{5.7824, 0.3907, 0.1422, 0.5362});
      desc.add("d1terms", std::vector{0.611, 1.039, 2.560});
      desc.add("d2terms", std::vector{2.103, 2.734, 2.835});
      desc.add("fs1Coefficients", std::vector{9.464, -9.054, -0.410});
      desc.add("fs2Coefficients", std::vector{-1.549, 1.985, -0.436});
      desc.add("fv1Coefficients", std::vector{-38.885, 425.007, -389.742});
      desc.add("fv2Coefficients", std::vector{-73.535, 83.211, -29.467});
      desc.add("LambdaSq", 9.733);
      desc.add("gamma", 2.148);
      return desc;
    }

  private:
    inline void eval() override {
      const double log1 = std::pow(log((lambda_sq_ + q2_) * inv_q20_), -gamma_);  // L(t=-q2) function in ref.

      // best fit parameterisation

      // scalar part
      auto Fs1{0.}, Fs2{0.};
      for (size_t i = 0; i < d1_terms_.size(); i++) {
        const auto d_1 = d1_terms_.at(i) + q2_;
        Fs1 += fs1_coefficients_.at(i) / d_1 * log1;
        Fs2 += fs2_coefficients_.at(i) / d_1 * log1;
      }

      // vector part
      auto d_2 = d2_terms_;
      for (auto& d : d_2)
        d += q2_;
      auto Fv1_term{0.}, Fv2_term{0.};
      for (size_t i = 0; i < d2_terms_.size(); i++) {
        Fv1_term += fv1_coefficients_.at(i) / d_2.at(i);
        Fv2_term += fv2_coefficients_.at(i) / d_2.at(i);
      }
      const auto log2 = std::pow(std::log((lambda_sq_ - 0.500) * inv_q20_), +gamma_),
                 log3 = std::pow(std::log((lambda_sq_ - 0.400) * inv_q20_), +gamma_);
      const auto Fv1 = (0.5 *
                            (rho1_coefficients_.at(0) * log2 +
                             rho1_coefficients_.at(1) * log3 * std::pow(1. + q2_ / rho1_coefficients_.at(2), -2)) /
                            (1. + q2_ / rho1_coefficients_.at(3)) +
                        Fv1_term) *
                       log1,
                 Fv2 = (0.5 *
                            (rho2_coefficients_.at(0) * log2 +
                             rho2_coefficients_.at(1) * log3 / (1. + q2_ / rho2_coefficients_.at(2))) /
                            (1. + q2_ / rho2_coefficients_.at(3)) +
                        Fv2_term) *
                       log1;

      const double F1 = Fv1 + Fs1, F2 = Fv2 + Fs2;
      setGEGM(F1 - tau(q2_) * F2, F1 + F2);
    }

    const std::vector<double> rho1_coefficients_, rho2_coefficients_;
    const std::vector<double> d1_terms_{}, d2_terms_{};
    const std::vector<double> fs1_coefficients_{}, fs2_coefficients_{};
    const std::vector<double> fv1_coefficients_{}, fv2_coefficients_{};
    const double inv_q20_;
    const double lambda_sq_;
    const double gamma_;
  };
}  // namespace cepgen::formfac
using cepgen::formfac::MergellEtAl;
REGISTER_FORMFACTORS("Mergell", MergellEtAl);
