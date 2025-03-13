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
#include "CepGen/Integration/FunctionalIntegrand.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"

using namespace cepgen;

FunctionalIntegrand::FunctionalIntegrand(const std::string& expr,
                                         const std::vector<std::string>& vars,
                                         const std::string& func_eval)
    : func_(FunctionalFactory::get().build(
          func_eval,
          ParametersList().set<std::string>("expression", expr).set<std::vector<std::string> >("variables", vars))) {
  CG_DEBUG("FunctionalIntegrand") << "Built a " << func_eval << " " << func_->variables().size()
                                  << "-dimensional functional with variables " << func_->variables() << " (" << vars
                                  << "): " << func_->expression() << ".";
}

double FunctionalIntegrand::eval(const std::vector<double>& x) {
  if (!func_)
    throw CG_FATAL("FunctionalIntegrand:eval") << "Functional object was not properly initialised!";
  if (x.size() != size())
    throw CG_FATAL("FunctionalIntegrand:eval")
        << "Invalid coordinates multiplicity: expected(" << size() << ") != received(" << x.size() << ")!";
  return (*func_)(x);
}

size_t FunctionalIntegrand::size() const {
  if (!func_)
    throw CG_FATAL("FunctionalIntegrand:eval") << "Functional object was not properly initialised!";
  return func_->variables().size();
}
