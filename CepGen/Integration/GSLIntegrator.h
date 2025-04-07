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

#ifndef CepGen_Integration_GSLIntegrator_h
#define CepGen_Integration_GSLIntegrator_h

#include "CepGen/Integration/Integrator.h"
#include "CepGen/Utils/GSLMonteFunctionWrapper.h"

namespace cepgen::utils {
  class RandomGenerator;
}  // namespace cepgen::utils

namespace cepgen {
  class GSLIntegrator : public Integrator {
  public:
    explicit GSLIntegrator(const ParametersList&);

    static ParametersDescription description();

  protected:
    void prepare(Integrand&, const std::vector<Limits>&);

    const std::unique_ptr<utils::RandomGenerator> random_generator_;  ///< GSL random number generator
    std::function<double(double*, size_t, void*)> function_;          ///< Functor wrapping GSL function footprint
    std::unique_ptr<gsl_monte_function> gsl_function_;                ///< Integrand, and its parameters
    std::vector<double> x_low_;                                       ///< Lower bounds to all integration variables
    std::vector<double> x_high_;                                      ///< Upper bounds to all integration variables
  };
}  // namespace cepgen

#endif
