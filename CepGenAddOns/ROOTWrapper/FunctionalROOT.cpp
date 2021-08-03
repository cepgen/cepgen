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

namespace cepgen {
  namespace utils {
    class FunctionalROOT final : public Functional {
    public:
      FunctionalROOT(const ParametersList& params) : Functional(params) {
        for (size_t i = 0; i < vars_.size(); ++i)
          func_.AddVariable(vars_[i], 0.);
        if (func_.Compile(expression_.c_str()) != 0)
          throw CG_ERROR("FunctionalROOT") << "Failed to define the function\n\t" << expression_;
      }
      static std::string description() { return "Plain old TFormula evaluator from ROOT"; }

      double eval(const std::vector<double>& x) const override {
        if (!func_.IsValid())
          throw CG_WARNING("FunctionalROOT") << "Cannot evaluate the invalid function at " << x << ".";
        return func_.EvalPar(values_.data());
      }

    private:
      TFormula func_;
    };
  }  // namespace utils
}  // namespace cepgen

REGISTER_FUNCTIONAL("ROOT", FunctionalROOT)
