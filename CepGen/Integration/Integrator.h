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

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  class Integrand;
  /// Monte-Carlo integration algorithm
  class Integrator : public NamedModule<std::string> {
  public:
    /// Integrator algorithm constructor
    Integrator(const ParametersList& params);

    static ParametersDescription parametersDescription();

    /// Specify the function to be integrated
    /// \param[in] integr Integrand object to be evaluated
    virtual void setIntegrand(Integrand& integr);
    /// Dimensional size of the phase space
    size_t size() const;

    /// Compute the function value at the given phase space point
    virtual double eval(const std::vector<double>& x) const;
    /// Generate a uniformly distributed (between 0 and 1) random number
    virtual double uniform() const;

    /// Perform the multidimensional Monte Carlo integration
    /// \param[out] result_ The cross section as integrated for the given phase space restrictions
    /// \param[out] abserr_ The uncertainty associated to the computed cross section
    virtual void integrate(double& result_, double& abserr_) = 0;

  protected:
    const unsigned long seed_;  ///< Random number generator seed
    int verbosity_;             ///< Integrator verbosity
    Integrand* integrand_;      ///< Integrand to be evaluated
    double result_{0.};         ///< Result of the last integration
    double err_result_{0.};     ///< Standard deviation for the last integration
    bool initialised_{false};   ///< Has the algorithm alreay been initialised?
    mutable std::default_random_engine rnd_gen_;
    mutable std::uniform_real_distribution<double> rnd_;
  };
}  // namespace cepgen

#endif
