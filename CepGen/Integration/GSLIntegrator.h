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

#ifndef CepGen_Integration_GSLIntegrator_h
#define CepGen_Integration_GSLIntegrator_h

#include <gsl/gsl_rng.h>

#include <functional>
#include <memory>

#include "CepGen/Integration/Integrator.h"
#include "CepGen/Utils/GSLFunctionsWrappers.h"

namespace cepgen {
  class GSLIntegrator : public Integrator {
  public:
    explicit GSLIntegrator(const ParametersList&);

    static ParametersDescription description();

    void setLimits(const std::vector<Limits>&) override;

  protected:
    void setIntegrand(Integrand&);
    /// A functor wrapping GSL's function footprint
    std::function<double(double*, size_t, void*)> funct_;
    /// GSL structure storing the function to be integrated by this
    /// integrator instance (along with its parameters)
    std::unique_ptr<gsl_monte_function> function_;
    std::vector<double> xlow_, xhigh_;
  };
}  // namespace cepgen

#endif
