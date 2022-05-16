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

#include <functional>

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  namespace utils {
    class Functional;
    class GSLIntegrator : public SteeredObject<GSLIntegrator> {
    public:
      explicit GSLIntegrator(const ParametersList& = ParametersList());

      static ParametersDescription description();

      //using Function1D = double (*)(double, void*);
      using Function1D = std::function<double(double)>;

      /// Evaluate the integral of a function at a given value
      /// \param[in] func function to integrate
      /// \param[in] xmin (optional) lower integration range
      /// \param[in] xmax (optional) upper integration range
      //double eval(const Functional& func) const;
      double eval(const Function1D& func, double xmin = INVALID, double xmax = INVALID) const;

    private:
      const Limits range_;
      static constexpr double INVALID = -999.999;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
