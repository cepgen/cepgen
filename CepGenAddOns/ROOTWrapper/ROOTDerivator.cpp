/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#include <TF1.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DerivatorFactory.h"
#include "CepGen/Utils/Derivator.h"
#include "CepGen/Utils/FunctionsWrappers.h"

namespace cepgen {
  class ROOTDerivator : public utils::Derivator {
  public:
    explicit ROOTDerivator(const ParametersList& params) : Derivator(params) {}

    static ParametersDescription description() {
      auto desc = Derivator::description();
      return desc;
    }

    /// Evaluate the derivative of a function at a given value
    /// \param[in] func function to derive
    /// \param[in] x coordinate
    /// \param[in] h (optional) step size ; if not provided, will use default algorithm value
    double derivate(const utils::Function1D& func, double x, double h = -1.) const override {
      auto rfunc = TF1(
          "cepgen_functional",
          [&func](double vars[1], double* pars) { return func(vars[0], static_cast<void*>(pars)); },
          0,
          1,
          0);
      return rfunc.Derivative(x, nullptr, h < 0. ? h_ : h);
    }

  private:
  };
}  // namespace cepgen

REGISTER_DERIVATOR("root", ROOTDerivator)
