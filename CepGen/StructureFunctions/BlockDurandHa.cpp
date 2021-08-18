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

#include <cassert>
#include <cmath>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    /// \f$F_2\f$ parameterisation from Block, Durand, and Ha \cite Block:2014kza
    class BlockDurandHa final : public Parameterisation {
    public:
      explicit BlockDurandHa(const ParametersList& params = ParametersList());
      BlockDurandHa& eval(double xbj, double q2) override;

    private:
      std::vector<double> a_, b_, c_;
      double n_;
      /// Effective mass spread parameter
      double lambda_;
      /// Asymptotic log-behaviour transition scale factor
      double mu2_;
      /// Squared effective mass (~VM mass)
      double m2_;
    };

    BlockDurandHa::BlockDurandHa(const ParametersList& params)
        : Parameterisation(params),
          a_(params.get<std::vector<double> >("a", {8.205e-4, -5.148e-2, -4.725e-3})),
          b_(params.get<std::vector<double> >("b", {2.217e-3, 1.244e-2, 5.958e-4})),
          c_(params.get<std::vector<double> >("c", {0.255e0, 1.475e-1})),
          n_(params.get<double>("n", 11.49)),
          lambda_(params.get<double>("lambda", 2.430)),
          mu2_(params.get<double>("mu2", 2.82)),
          m2_(params.get<double>("m2", 0.753)) {
      assert(a_.size() == 3);
      assert(b_.size() == 3);
      assert(c_.size() == 2);
    }

    BlockDurandHa& BlockDurandHa::eval(double xbj, double q2) {
      if (q2 <= 0) {
        F2 = 0.;
        return *this;
      }

      const double tau = q2 / (q2 + mu2_);
      const double xl = log1p(q2 / mu2_);
      const double xlx = log(tau / xbj);

      const double A = a_[0] + a_[1] * xl + a_[2] * xl * xl;
      const double B = b_[0] + b_[1] * xl + b_[2] * xl * xl;
      const double C = c_[0] + c_[1] * xl;
      const double D = q2 * (q2 + lambda_ * m2_) / pow(q2 + m2_, 2);

      F2 = D * pow(1. - xbj, n_) * (C + A * xlx + B * xlx * xlx);

      return *this;
    }
  }  // namespace strfun
}  // namespace cepgen

REGISTER_STRFUN(BlockDurandHa, strfun::BlockDurandHa)
