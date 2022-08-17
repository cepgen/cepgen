/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#ifndef CepGen_Utils_GSLIntegrator_h
#define CepGen_Utils_GSLIntegrator_h

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Utils/GSLFunctionsWrappers.h"

namespace cepgen {
  namespace utils {
    class GSLIntegrator : public SteeredObject<GSLIntegrator> {
    public:
      explicit GSLIntegrator(const ParametersList& = ParametersList());

      static ParametersDescription description();

      /// Evaluate the derivative of a function at a given value
      /// \param[in] func function to derive
      /// \param[in] x coordinate
      /// \param[in] h (optional) step size ; if not provided, will use default algorithm value
      double eval(const std::function<double(double)>& func, double xmin = INVALID, double xmax = INVALID) const {
        return eval(Function1D(func), xmin, xmax);
      }
      /// Evaluate the integral of a function at a given value
      /// \param[in] func function to integrate
      /// \param[in] xmin (optional) lower integration range
      /// \param[in] xmax (optional) upper integration range
      double eval(const Function1D& func, double xmin = INVALID, double xmax = INVALID) const;
      /// Evaluate the integral of a function at a given value
      /// \param[in] func function to integrate
      /// \param[in] obj parameters object
      /// \param[in] xmin (optional) lower integration range
      /// \param[in] xmax (optional) upper integration range
      double eval(const Function1D& func, void* obj, double xmin = INVALID, double xmax = INVALID) const;

    private:
      enum struct Mode { Fixed = 0, QNG = 1, QAG = 2, QAGS = 3, QAWC = 4 };
      enum struct FixedType {
        Legendre = 0,
        Chebyshev = 1,
        Gegenbauer = 2,
        Jacobi = 3,
        Laguerre = 4,
        Hermite = 5,
        Exponential = 6,
        Rational = 7,
        Chebyshev2 = 8
      };
      double eval(const gsl_function*, double xmin = INVALID, double xmax = INVALID) const;
      const Limits range_;
      const Mode mode_;
      const FixedType fixed_type_;
      const int nodes_;
      const double alpha_, beta_;
      const size_t limit_;
      const double epsabs_, epsrel_;
      const ParametersList func_params_;
      static constexpr double INVALID = -999.999;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
