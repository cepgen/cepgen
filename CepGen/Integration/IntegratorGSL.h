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

#ifndef CepGen_Integration_IntegratorGSL_h
#define CepGen_Integration_IntegratorGSL_h

#include <gsl/gsl_rng.h>

#include <functional>
#include <memory>

#include "CepGen/Integration/Integrator.h"
#include "CepGen/Utils/GSLFunctionsWrappers.h"

namespace cepgen {
  class IntegratorGSL : public Integrator {
  public:
    explicit IntegratorGSL(const ParametersList&);

    static ParametersDescription description();

    double uniform(const Limits&) const override;
    void setLimits(const std::vector<Limits>&) override;

  protected:
    void setIntegrand(Integrand&);
    /// A functor wrapping GSL's function footprint
    std::function<double(double*, size_t, void*)> funct_;
    /// GSL structure storing the function to be integrated by this
    /// integrator instance (along with its parameters)
    std::unique_ptr<gsl_monte_function> function_;
    /// A deleter object for GSL's random number generator
    struct gsl_rng_deleter {
      /// Destructor method for the random number generator service
      inline void operator()(gsl_rng* rng) { gsl_rng_free(rng); }
    };
    /// Instance of random number generator service
    std::unique_ptr<gsl_rng, gsl_rng_deleter> gsl_rng_;
    std::vector<double> xlow_, xhigh_;
  };
}  // namespace cepgen

#endif
