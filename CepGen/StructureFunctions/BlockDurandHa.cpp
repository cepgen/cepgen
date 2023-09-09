/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    /// \f$F_2\f$ parameterisation from Block, Durand, and Ha \cite Block:2014kza
    class BlockDurandHa final : public Parameterisation {
    public:
      explicit BlockDurandHa(const ParametersList& params)
          : Parameterisation(params),
            a_(steer<std::vector<double> >("a")),
            b_(steer<std::vector<double> >("b")),
            c_(steer<std::vector<double> >("c")),
            n_(steer<double>("n")),
            lambda_(steer<double>("lambda")),
            mu2_(steer<double>("mu2")),
            m2_(steer<double>("m2")) {
        if (a_.size() != 3)
          throw CG_FATAL("BlockDurandHa") << "Parameter 'a' should have 3 components! Parsed " << a_ << ".";
        if (b_.size() != 3)
          throw CG_FATAL("BlockDurandHa") << "Parameter 'b' should have 3 components! Parsed " << b_ << ".";
        if (c_.size() != 2)
          throw CG_FATAL("BlockDurandHa") << "Parameter 'c' should have 3 components! Parsed " << c_ << ".";
      }

      void eval() override {
        const double tau = args_.q2 / (args_.q2 + mu2_);
        const double xl = log1p(args_.q2 / mu2_);
        const double xlx = log(tau / args_.xbj);

        const double A = a_.at(0) + a_.at(1) * xl + a_.at(2) * xl * xl;
        const double B = b_.at(0) + b_.at(1) * xl + b_.at(2) * xl * xl;
        const double C = c_.at(0) + c_.at(1) * xl;
        const double D = args_.q2 * (args_.q2 + lambda_ * m2_) / pow(args_.q2 + m2_, 2);

        setF2(D * pow(1. - args_.xbj, n_) * (C + A * xlx + B * xlx * xlx));
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Block-Durand-Ha (continuum)");
        desc.add<std::vector<double> >("a", {8.205e-4, -5.148e-2, -4.725e-3});
        desc.add<std::vector<double> >("b", {2.217e-3, 1.244e-2, 5.958e-4});
        desc.add<std::vector<double> >("c", {0.255e0, 1.475e-1});
        desc.add<double>("n", 11.49);
        desc.add<double>("lambda", 2.430);
        desc.add<double>("mu2", 2.82);
        desc.add<double>("m2", 0.753);
        return desc;
      }

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
  }  // namespace strfun
}  // namespace cepgen
using cepgen::strfun::BlockDurandHa;
REGISTER_STRFUN(13, BlockDurandHa);
