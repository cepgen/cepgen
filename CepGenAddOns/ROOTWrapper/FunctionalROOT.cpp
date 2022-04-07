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

#include <TFormula.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    class FunctionalROOT final : public Functional {
    public:
      explicit FunctionalROOT(const ParametersList&);
      double eval() const override;

      static ParametersDescription description();

    private:
      TFormula func_;
    };

    FunctionalROOT::FunctionalROOT(const ParametersList& params) : Functional(params) {
      for (auto& var : vars_)
        func_.AddVariable(var, 0.);
      auto expr = expression_;
      expr = utils::replace_all(expr, {{"min(", "TMath::Min("}, {"max(", "TMath::Max("}});
      if (func_.Compile(expr.c_str()) != 0)
        throw CG_ERROR("FunctionalROOT") << "Failed to define the function\n\t" << expression_;
      CG_DEBUG("FunctionalROOT") << "Successfully defined a dimension-" << vars_.size() << " function with arguments "
                                 << vars_ << ": " << expr << ".";
    }

    double FunctionalROOT::eval() const {
      if (!func_.IsValid())
        throw CG_WARNING("FunctionalROOT") << "Cannot evaluate the invalid function at " << values_ << ".";
      return func_.EvalPar(values_.data());
    }

    ParametersDescription FunctionalROOT::description() {
      auto desc = Functional::description();
      desc.setDescription("Plain old TFormula evaluator from ROOT");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_FUNCTIONAL("ROOT", FunctionalROOT)
