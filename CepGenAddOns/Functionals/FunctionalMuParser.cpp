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

#include <muParser.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"

namespace cepgen {
  namespace utils {
    class FunctionalMuParser final : public Functional {
    public:
      explicit FunctionalMuParser(const ParametersList&);

      static ParametersDescription parametersDescription();

      double eval(const std::vector<double>&) const override;

    private:
      mu::Parser parser_;
    };

    FunctionalMuParser::FunctionalMuParser(const ParametersList& params) : Functional(params) {
      try {
        for (size_t i = 0; i < vars_.size(); ++i)
          parser_.DefineVar(vars_[i], &values_[i]);
        parser_.SetExpr(expression_);
      } catch (const mu::Parser::exception_type& e) {
        throw CG_ERROR("FunctionalMuParser")
            << "Failed to define the function\n\t" << expression_ << "\n\t" << std::string(e.GetPos(), '-') + "^"
            << "\n\t" << e.GetMsg();
      }
    }

    double FunctionalMuParser::eval(const std::vector<double>& x) const {
      try {
        values_ = x;
        return parser_.Eval();
      } catch (const mu::Parser::exception_type& e) {
        throw CG_WARNING("FunctionalMuParser")
            << "Failed to evaluate the function\n\t" << expression_ << "\n\t" << std::string(e.GetPos(), '-') + "^"
            << "\n\t" << e.GetMsg();
      }
    }

    ParametersDescription FunctionalMuParser::parametersDescription() {
      auto desc = Functional::parametersDescription();
      desc.setDescription("MuParser functional evaluator");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_FUNCTIONAL("MuParser", FunctionalMuParser)
