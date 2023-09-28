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

#ifndef CepGen_Integration_Integrator_h
#define CepGen_Integration_Integrator_h

#include <random>
#include <vector>

#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Utils/Value.h"

namespace cepgen {
  class Integrand;
  /// Monte-Carlo integration algorithm
  class Integrator : public NamedModule<std::string> {
  public:
    /// Integrator algorithm constructor
    explicit Integrator(const ParametersList& params);

    static ParametersDescription description();

    /// Ensure the integration bounds are properly set
    void checkLimits(const Integrand&);
    /// Specify the variables limits on integration
    virtual void setLimits(const std::vector<Limits>& limits) { limits_ = limits; }

    /// Compute the function value at the given phase space point
    virtual double eval(Integrand&, const std::vector<double>&) const;
    /// Generate a uniformly distributed (between 0 and 1) random number
    virtual double uniform(double min = 0., double max = 1.) const;

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
    const unsigned long seed_;    ///< Random number generator seed
    int verbosity_;               ///< Integrator verbosity
    std::vector<Limits> limits_;  ///< List of per-variable integration limits
    mutable std::default_random_engine rnd_gen_;
    mutable std::uniform_real_distribution<double> rnd_;
  };
}  // namespace cepgen

#endif
