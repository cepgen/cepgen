/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

#ifndef CepGen_Utils_Derivator_h
#define CepGen_Utils_Derivator_h

#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Utils/FunctionWrapper.h"

namespace cepgen::utils {
  class Derivator : public NamedModule<Derivator> {
  public:
    explicit Derivator(const ParametersList& params) : NamedModule(params), h_(steer<double>("h")) {}

    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.add<double>("h", 1.e-2).setDescription("step size");
      return desc;
    }

    /// Evaluate the derivative of a function at a given value
    /// \param[in] func function to derive
    /// \param[in] x coordinate
    /// \param[in] h (optional) step size ; if not provided, will use default algorithm value
    inline double derivate(const std::function<double(double)>& func, double x, double h = -1.) const {
      return derivate(FunctionWrapper(func), x, h);
    }
    /// Evaluate the derivative of a function at a given value
    /// \param[in] func function to derive
    /// \param[in] x coordinate
    /// \param[in] h (optional) step size ; if not provided, will use default algorithm value
    virtual double derivate(const FunctionWrapper& func, double x, double h = -1.) const = 0;

  protected:
    const double h_;
  };
}  // namespace cepgen::utils

#endif
