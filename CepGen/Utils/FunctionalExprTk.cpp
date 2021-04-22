#include "CepGen/Utils/Functional.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Core/Exception.h"

#include <exprtk.hpp>

namespace cepgen {
  namespace utils {
    class FunctionalExprTk : public Functional {
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
      parser_.compile(expression_, expr_);
    }

    double FunctionalExprTk::eval(const std::vector<double>& x) const { return expr_.value(); }
  }  // namespace utils
}  // namespace cepgen

REGISTER_FUNCTIONAL("ExprTk", FunctionalExprTk)
