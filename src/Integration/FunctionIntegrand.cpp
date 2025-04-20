/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

using namespace cepgen;

FunctionIntegrand::FunctionIntegrand(size_t num_dimensions,
                                     const std::function<double(const std::vector<double>&)>& function)
    : function_(function), num_dimensions_(num_dimensions) {}

double FunctionIntegrand::eval(const std::vector<double>& coordinates) {
  if (coordinates.size() != size())
    throw CG_FATAL("FunctionIntegrand:eval")
        << "Invalid coordinates multiplicity: expected(" << size() << ") != received(" << coordinates.size() << ")!";
  const auto weight = function_(coordinates);  // calculate weight for the phase space point to probe
  CG_DEBUG_LOOP("FunctionIntegrand:eval")
      << "f value for dim-" << coordinates.size() << " point " << coordinates << ": " << weight << ".";
  return weight;
}
