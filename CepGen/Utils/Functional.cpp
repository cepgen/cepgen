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
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    Functional::Functional(const ParametersList& params)
        : NamedModule(params),
          vars_orig_(steer<std::vector<std::string> >("variables")),
          expression_orig_(steer<std::string>("expression")),
          vars_(vars_orig_),
          expression_(expression_orig_),
          values_(vars_.size()) {
      for (size_t i = 0; i < vars_.size(); ++i) {
        replace_all(vars_.at(i), "(", "_");
        replace_all(vars_.at(i), ")", "_");
        replace_all(expression_, vars_orig_.at(i), vars_.at(i));
      }
    }

    double Functional::operator()(double x) const {
      if (vars_orig_.size() != 1)
        throw CG_FATAL("Functional") << "This function only works with single-dimensional functions!";
      return operator()(std::vector<double>{x});
    }

    double Functional::operator()(const std::vector<double>& x) const {
      values_ = x;
      return eval(x);
    }

    ParametersDescription Functional::description() {
      auto desc = ParametersDescription();
      desc.setDescription("Unnamed functional evaluator");
      desc.add<std::vector<std::string> >("variables", {}).setDescription("List of variables to evaluate");
      desc.add<std::string>("expression", "").setDescription("Functional expression");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen
