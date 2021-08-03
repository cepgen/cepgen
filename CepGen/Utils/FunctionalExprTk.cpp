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

#include <exprtk.hpp>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    class FunctionalExprTk final : public Functional {
    public:
      explicit FunctionalExprTk(const ParametersList&);
      double eval(const std::vector<double>&) const override;

    private:
      exprtk::symbol_table<double> symbols_;
      exprtk::expression<double> expr_;
      exprtk::parser<double> parser_;
    };

    FunctionalExprTk::FunctionalExprTk(const ParametersList& params) : Functional(params) {
      for (size_t i = 0; i < vars_.size(); ++i)
        symbols_.add_variable(vars_[i], values_[i]);
      symbols_.add_constants();
      expr_.register_symbol_table(symbols_);
      replace_all(expression_, "**", "^");
      if (!parser_.compile(expression_, expr_))
        throw CG_WARNING("FunctionalExprTk") << "Failed to compile expression \"" << expression() << "\".";
    }

    double FunctionalExprTk::eval(const std::vector<double>&) const { return expr_.value(); }
  }  // namespace utils
}  // namespace cepgen

REGISTER_FUNCTIONAL("ExprTk", FunctionalExprTk)
