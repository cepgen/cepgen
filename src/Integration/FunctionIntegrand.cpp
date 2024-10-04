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

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/FunctionIntegrand.h"

namespace cepgen {
  FunctionIntegrand::FunctionIntegrand(size_t num_dimensions,
                                       const std::function<double(const std::vector<double>&)>& func)
      : function_(func), num_dimensions_(num_dimensions) {}

  double FunctionIntegrand::eval(const std::vector<double>& x) {
    if (x.size() != size())
      throw CG_FATAL("FunctionIntegrand:eval")
          << "Invalid coordinates multiplicity: expected(" << size() << ") != received(" << x.size() << ")!";

    //--- calculate weight for the phase space point to probe
    double weight = function_(x);

    //--- a bit of useful debugging
    CG_DEBUG_LOOP("FunctionIntegrand:eval")
        << "f value for dim-" << x.size() << " point " << x << ": " << weight << ".";

    return weight;
  }
}  // namespace cepgen
