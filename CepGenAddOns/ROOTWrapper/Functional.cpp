/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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
  namespace root {
    /// Functional evaluator defined from a ROOT TFormula
    class Functional final : public cepgen::utils::Functional {
    public:
      explicit Functional(const ParametersList& params) : cepgen::utils::Functional(params) {
        for (auto& var : vars_)
          func_.AddVariable(var, 0.);
        auto expr = expression_;
        expr = utils::replaceAll(expr, {{"min(", "TMath::Min("}, {"max(", "TMath::Max("}});
        if (func_.Compile(expr.c_str()) != 0)
          throw CG_ERROR("root:Functional") << "Failed to define the function\n\t" << expression_;
        CG_DEBUG("root:Functional") << "Successfully defined a dimension-" << vars_.size()
                                    << " function with arguments " << vars_ << ": " << expr << ".";
      }
      inline double eval() const override {
        if (!func_.IsValid())
          throw CG_WARNING("root:Functional") << "Cannot evaluate the invalid function at " << values_ << ".";
        return func_.EvalPar(values_.data());
      }

      inline static ParametersDescription description() {
        auto desc = cepgen::utils::Functional::description();
        desc.setDescription("Plain old TFormula evaluator from ROOT");
        return desc;
      }

    private:
      TFormula func_;
    };
  }  // namespace root
}  // namespace cepgen
using cepgen::root::Functional;
REGISTER_FUNCTIONAL("root", Functional);
