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

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/TimeKeeper.h"

namespace cepgen {
  Integrand::~Integrand() { CG_DEBUG("Integrand") << "Destructor called"; }

  double FunctionIntegrand::eval(const std::vector<double>& x) {
    if (x.size() != size())
      throw CG_FATAL("Integrand:eval") << "Invalid coordinates multiplicity: expected(" << size() << ") != received("
                                       << x.size() << ")!";

    //--- calculate weight for the phase space point to probe
    double weight = function_(x);

    //--- a bit of useful debugging
    CG_DEBUG_LOOP("Integrand") << "f value for dim-" << x.size() << " point " << x << ": " << weight << ".";

    return weight;
  }
}  // namespace cepgen
