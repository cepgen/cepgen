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
  Integrand::Integrand(const Parameters* params) : params_(params), tmr_(new utils::Timer) {}

  Integrand::~Integrand() { CG_DEBUG("Integrand") << "Destructor called"; }

  void Integrand::setFunction(size_t ndim, const std::function<double(const std::vector<double>&)>& func) {
    gen_integr_.ndim = ndim;
    gen_integr_.function = func;
  }

  double Integrand::eval(const std::vector<double>& x) {
    CG_TICKER(const_cast<Parameters*>(params_)->timeKeeper());

    if (x.size() != size())
      throw CG_FATAL("Integrand:eval") << "Invalid coordinates multiplicity: expected(" << size() << ") != received("
                                       << x.size() << ")!";

    //--- start the timer
    tmr_->reset();

    //--- specify the phase space point to probe and calculate weight
    double weight = gen_integr_.function(x);

    //--- a bit of useful debugging
    CG_DEBUG_LOOP("Integrand") << "f value for dim-" << x.size() << " point " << x << ": " << weight << ".";

    return weight;
  }
}  // namespace cepgen
