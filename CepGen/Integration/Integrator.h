/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#ifndef CepGen_Integration_Integrator_h
#define CepGen_Integration_Integrator_h

#include <functional>
#include <memory>

#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/Value.h"

namespace cepgen {
  class Integrand;
  /// Monte-Carlo integration algorithm
  class Integrator : public NamedModule<Integrator> {
  public:
    /// Integrator algorithm constructor
    explicit Integrator(const ParametersList& params);

    static ParametersDescription description();

    void checkLimits(const Integrand&);  ///< Ensure the integration bounds are properly set
    /// Set variables integration limits
    virtual void setLimits(const std::vector<Limits>& limits) { limits_ = limits; }

    virtual double eval(Integrand&, const std::vector<double>&) const;  ///< Compute function value at one point

    /// Generate a uniformly distributed (between 0 and 1) random number
    virtual double uniform(const Limits& = {0., 1.}) const;

    /// Perform the multidimensional Monte Carlo integration
    /// \param[out] result integral computed over the full phase space
    virtual Value integrate(Integrand& result) = 0;
    /// Perform an integration with a given functional and a given set of parameters
    static Value integrate(const std::function<double(const std::vector<double>&)>&, const ParametersList&, size_t);
    /// Perform an integration with a given functional and a given set of parameters
    static Value integrate(const std::function<double(const std::vector<double>&)>&,
                           const ParametersList&,
                           const std::vector<Limits>&);

  protected:
    const std::unique_ptr<utils::RandomGenerator> random_number_generator_;
    int verbosity_;               ///< Integrator verbosity
    std::vector<Limits> limits_;  ///< List of per-variable integration limits
  };
}  // namespace cepgen

#endif
