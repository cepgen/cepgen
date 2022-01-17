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

#ifndef CepGen_Utils_GSLDerivator_h
#define CepGen_Utils_GSLDerivator_h

#include <functional>

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  namespace utils {
    class GSLDerivator : public SteeredObject<GSLDerivator> {
    public:
      explicit GSLDerivator(const ParametersList&);

      static ParametersDescription description();

      /// A one-dimensional function to evaluate
      typedef std::function<double(double)> Function1D;

      /// Evaluate the derivative of a function at a given value
      /// \param[in] x coordinate
      /// \param[in] h (optional) step size ; if not provided, will use default algorithm value
      double eval(const Function1D& func, double x, double h = -1.) const;

      enum struct Mode {
        central,  ///< adaptive central difference algorithm
        forward,  ///< adaptive forward difference algorithm
        backward  ///< adaptive backward difference algorithm
      };

    private:
      const Mode mode_;
      const double h_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
