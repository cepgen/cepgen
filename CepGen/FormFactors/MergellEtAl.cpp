/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"

namespace cepgen {
  namespace formfac {
    /// \cite Mergell:1995bf
    class MergellEtAl final : public Parameterisation {
    public:
      explicit MergellEtAl(const ParametersList&);

      static ParametersDescription description();

    private:
      FormFactors compute(double q2) override;

      const double a1rho_, a2rho_, b1rho_, b2rho_, c1rho_, c2rho_, d1rho_, d2rho_;
      const double inv_q20_;
      const double lambda_sq_;
      const double gamma_;
    };

    MergellEtAl::MergellEtAl(const ParametersList& params)
        : Parameterisation(params),
          a1rho_(steer<double>("a1rho")),
          a2rho_(steer<double>("a2rho")),
          b1rho_(steer<double>("b1rho")),
          b2rho_(steer<double>("b2rho")),
          c1rho_(steer<double>("c1rho")),
          c2rho_(steer<double>("c2rho")),
          d1rho_(steer<double>("d1rho")),
          d2rho_(steer<double>("d2rho")),
          inv_q20_(steer<double>("q20inv")),
          lambda_sq_(steer<double>("lambdaSq")),
          gamma_(steer<double>("gamma")) {}

    FormFactors MergellEtAl::compute(double q2) {
      const double log1 = std::pow(log((lambda_sq_ + q2) * inv_q20_), -gamma_);  // L(t=-q2) function in ref.

      // best fit parameterisation
      const double d1_1 = 0.611 + q2, d2_1 = 1.039 + q2, d3_1 = 2.560 + q2;
      const double Fs1 = (9.464 / d1_1 - 9.054 / d2_1 - 0.410 / d3_1) * log1;
      const double Fs2 = (-1.549 / d1_1 + 1.985 / d2_1 - 0.436 / d3_1) * log1;

      const double log2 = std::pow(log((lambda_sq_ - 0.500) * inv_q20_), +gamma_);
      const double log3 = std::pow(log((lambda_sq_ - 0.400) * inv_q20_), +gamma_);

      const double d1_2 = 2.103 + q2, d2_2 = 2.734 + q2, d3_2 = 2.835 + q2;
      const double Fv1 = (0.5 * (a1rho_ * log2 + b1rho_ * log3 * std::pow(1. + q2 / c1rho_, -2)) / (1. + q2 / d1rho_) -
                          38.885 / d1_2 + 425.007 / d2_2 - 389.742 / d3_2) *
                         log1;
      const double Fv2 = (0.5 * (a2rho_ * log2 + b2rho_ * log3 / (1. + q2 / c2rho_)) / (1. + q2 / d2rho_) -
                          73.535 / d1_2 + 83.211 / d2_2 - 29.467 / d3_2) *
                         log1;

      const double F1 = Fv1 + Fs1, F2 = Fv2 + Fs2;

      FormFactors out;
      out.GE = F1 - tau(q2) * F2;
      out.GM = F1 + F2;
      return out;
    }

    ParametersDescription MergellEtAl::description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Mergell et al.");
      desc.add<double>("a1rho", 1.0317);
      desc.add<double>("a2rho", 5.7824);
      desc.add<double>("b1rho", 0.0875);
      desc.add<double>("b2rho", 0.3907);
      desc.add<double>("c1rho", 0.3176);
      desc.add<double>("c2rho", 0.1422);
      desc.add<double>("d1rho", 0.5496);
      desc.add<double>("d2rho", 0.5362);
      desc.add<double>("q20inv", 1. / 0.35);
      desc.add<double>("lambdaSq", 9.733);
      desc.add<double>("gamma", 2.148);
      return desc;
    }
  }  // namespace formfac
}  // namespace cepgen

REGISTER_FORMFACTORS("Mergell", MergellEtAl)
