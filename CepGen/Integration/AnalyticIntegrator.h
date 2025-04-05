/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#ifndef CepGen_Integration_AnalyticIntegrator_h
#define CepGen_Integration_AnalyticIntegrator_h

#include <functional>

#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen::utils {
  class FunctionWrapper;
}

namespace cepgen {
  /// Analytic (functional) integration algorithm
  class AnalyticIntegrator : public NamedModule<AnalyticIntegrator> {
  public:
    explicit AnalyticIntegrator(const ParametersList&);

    static ParametersDescription description();

    /// Evaluate the integral of a function at a given value
    /// \param[in] integrand function to integrate
    /// \param[in] range (optional) integration range
    double integrate(const std::function<double(double)>& integrand, const Limits& range = {}) const;
    /// Evaluate the integral of a function at a given value
    /// \param[in] integrand function to integrate
    /// \param[in] parameters specific parameters object
    /// \param[in] range (optional) integration range
    template <typename T>
    double integrate(const utils::FunctionWrapper& integrand, const T& parameters, const Limits& range = {}) const {
      return run(integrand, const_cast<void*>(reinterpret_cast<const void*>(&parameters)), range);
    }

  protected:
    /// Evaluate the integral of a function at a given value
    /// \param[in] integrand function to integrate
    /// \param[in] parameters (optional) parameters object
    /// \param[in] range (optional) integration range
    virtual double run(const utils::FunctionWrapper& integrand,
                       void* parameters = nullptr,
                       const Limits& range = {}) const = 0;

    const Limits range_;
    const ParametersList integrand_parameters_;
    const int verbosity_;  ///< Integrator verbosity
  };
}  // namespace cepgen

#endif
