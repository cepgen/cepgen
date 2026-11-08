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

FunctionalIntegrand::FunctionalIntegrand(const std::string& expression,
                                         const std::vector<std::string>& variables,
                                         const std::string& functional_evaluator)
    : functional_(FunctionalFactory::get().build(
          functional_evaluator, ParametersList().set("expression", expression).set("variables", variables))) {
  CG_DEBUG("FunctionalIntegrand") << "Built a " << functional_evaluator << " " << functional_->variables().size()
                                  << "-dimensional functional with variables " << functional_->variables() << " ("
                                  << variables << "): " << functional_->expression() << ".";
}

double FunctionalIntegrand::eval(const std::vector<double>& coordinates) {
  if (!functional_)
    throw CG_FATAL("FunctionalIntegrand:eval") << "Functional object was not properly initialised!";
  if (coordinates.size() != size())
    throw CG_FATAL("FunctionalIntegrand:eval")
        << "Invalid coordinates multiplicity: expected(" << size() << ") != received(" << coordinates.size() << ")!";
  return (*functional_)(coordinates);
}

size_t FunctionalIntegrand::size() const {
  if (!functional_)
    throw CG_FATAL("FunctionalIntegrand:eval") << "Functional object was not properly initialised!";
  return functional_->variables().size();
}
