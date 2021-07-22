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
