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

namespace cepgen {
  class Integrand;
  /// Monte-Carlo integration algorithm
  class Integrator : public NamedModule<std::string> {
  public:
    /// Integrator algorithm constructor
    Integrator(const ParametersList& params);

    static ParametersDescription description();

    /// Specify the function to be integrated
    /// \param[in] integr Integrand object to be evaluated
    virtual void setIntegrand(Integrand& integr);
    /// Specify the variables limits on integration
    void setLimits(const std::vector<Limits>& limits) { limits_ = limits; }
    /// Dimensional size of the phase space
    size_t size() const;

    /// Compute the function value at the given phase space point
    virtual double eval(const std::vector<double>& x) const;
    /// Generate a uniformly distributed (between 0 and 1) random number
    virtual double uniform(double min = 0., double max = 1.) const;

    /// Perform the multidimensional Monte Carlo integration
    /// \param[out] result integral computed over the full phase space
    /// \param[out] abserr uncertainty associated to the computed integral
    virtual void integrate(double& result, double& abserr) = 0;
    /// Perform an integration with no use of the numerical error
    /// \return the integral computed over the full phase space
    double integrate();
    /// Perform an integration with a given functional and a given set of parameters
    static double integrate(const std::function<double(const std::vector<double>&)>&, const ParametersList&, size_t);
    /// Perform an integration with a given functional and a given set of parameters
    static double integrate(const std::function<double(const std::vector<double>&)>&,
                            const ParametersList&,
                            const std::vector<Limits>&);

  protected:
    const unsigned long seed_;       ///< Random number generator seed
    int verbosity_;                  ///< Integrator verbosity
    Integrand* integrand_{nullptr};  ///< Integrand to be evaluated
    std::vector<Limits> limits_;     ///< List of per-variable integration limits
    double result_{0.};              ///< Result of the last integration
    double err_result_{0.};          ///< Standard deviation for the last integration
    bool initialised_{false};        ///< Has the algorithm alreay been initialised?
    mutable std::default_random_engine rnd_gen_;
    mutable std::uniform_real_distribution<double> rnd_;
  };
}  // namespace cepgen

#endif
