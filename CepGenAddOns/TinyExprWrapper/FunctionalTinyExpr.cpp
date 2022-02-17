/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include <tinyexpr.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    class FunctionalTinyExpr final : public Functional {
    public:
      explicit FunctionalTinyExpr(const ParametersList&);
      double eval() const override;

      static ParametersDescription description();

    private:
      struct te_deleter {
        void operator()(te_expr* expr) { te_free(expr); }
      };
      std::unique_ptr<te_expr, te_deleter> eval_;
    };

    FunctionalTinyExpr::FunctionalTinyExpr(const ParametersList& params) : Functional(params) {
      std::vector<te_variable> te_vars;
      for (size_t i = 0; i < vars_.size(); ++i)
        te_vars.emplace_back(te_variable{vars_.at(i).c_str(), &values_.at(i), TE_VARIABLE, nullptr});
      auto expr = expression_;
      expr = utils::replace_all(expr, {{"**", "^"}});
      int error;
      eval_.reset(te_compile(expr.c_str(), te_vars.data(), vars_.size(), &error));
      if (!eval_) {
        const std::string pre_syntax_err = "A syntax error was detected in the expression \"";
        const std::string postfix = (expr != expression_) ? " (adapted from \"" + expression_ + "\")" : "";
        throw CG_ERROR("FunctionalTinyExpr") << "Evaluator was not properly initialised.\n"
                                             << pre_syntax_err << expr << "\"" << postfix << "\n"
                                             << std::string(pre_syntax_err.size() + error - 1, ' ') + "^";
      }
    }

    double FunctionalTinyExpr::eval() const { return te_eval(eval_.get()); }

    ParametersDescription FunctionalTinyExpr::description() {
      auto desc = Functional::description();
      desc.setDescription("TinyExpr evaluator");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_FUNCTIONAL("tinyexpr", FunctionalTinyExpr)
