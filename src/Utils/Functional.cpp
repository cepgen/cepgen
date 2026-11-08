/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;
using namespace cepgen::utils;
using namespace std::string_literals;

Functional::Functional(const ParametersList& params)
    : NamedModule(params),
      vars_orig_(steer<std::vector<std::string> >("variables")),
      expression_orig_(steer<std::string>("expression")),
      vars_(vars_orig_),
      expression_(expression_orig_),
      values_(vars_.size()) {
  for (size_t i = 0; i < vars_.size(); ++i) {
    vars_.at(i) = sanitise(vars_.at(i));
    replaceAll(expression_, vars_orig_.at(i), vars_.at(i));
  }
}

double Functional::operator()(double x) const {
  if (vars_.size() != 1)
    throw CG_FATAL("Functional") << "This function only works with single-dimensional functions!";
  values_[0] = x;
  return eval();
}

double Functional::operator()(const std::vector<double>& x) const {
  if (vars_.size() != x.size())
    throw CG_FATAL("Functional") << "Invalid number of variables fed to the evaluator! Expecting " << vars_.size()
                                 << ", got " << x.size() << ".";
  values_ = x;
  return eval();
}

ParametersList Functional::fromExpression(const std::string& expr, const std::vector<std::string>& vars) {
  return ParametersList().set("expression", expr).set("variables", vars);
}

ParametersDescription Functional::description() {
  auto desc = ParametersDescription();
  desc.setDescription("Unnamed functional evaluator");
  desc.add("variables", std::vector<std::string>{}).setDescription("List of variables to evaluate");
  desc.add("expression", ""s).setDescription("Functional expression");
  return desc;
}
