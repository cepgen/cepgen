/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025  Laurent Forthomme
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

#include <fparser.hh>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

namespace cepgen::fparser {
  class Functional final : public cepgen::utils::Functional {
  public:
    explicit Functional(const ParametersList& params)
        : cepgen::utils::Functional(params), function_parser_(new FunctionParser) {
      function_parser_->Parse(expression_, utils::merge(vars_, ","));
      //throw CG_ERROR("atmsp:Functional") << "Evaluator was not properly initialised. ATMSP error:\n" << parser.errMessage(err);
    }
    double eval() const override { return function_parser_->Eval(values_.data()); }

    static ParametersDescription description() {
      auto desc = cepgen::utils::Functional::description();
      desc.setDescription("fparser evaluator");
      return desc;
    }

  private:
    const std::unique_ptr<FunctionParser> function_parser_;
  };
}  // namespace cepgen::fparser
using FunctionParserFunctional = cepgen::fparser::Functional;
REGISTER_FUNCTIONAL("fparser", FunctionParserFunctional);
