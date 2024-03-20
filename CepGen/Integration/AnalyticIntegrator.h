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

#ifndef CepGen_Integration_AnalyticIntegrator_h
#define CepGen_Integration_AnalyticIntegrator_h

#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  namespace utils {
    class Function1D;
  }
  /// Analytic (functional) integration algorithm
  class AnalyticIntegrator : public NamedModule<AnalyticIntegrator, std::string> {
  public:
    explicit AnalyticIntegrator(const ParametersList& params);  ///< Integrator algorithm constructor

    static ParametersDescription description();

    /// Evaluate the integral of a function at a given value
    /// \param[in] func function to integrate
    /// \param[in] range (optional) integration range
    double integrate(const std::function<double(double)>& func, const Limits& range = {}) const;
    /// Evaluate the integral of a function at a given value
    /// \param[in] func function to integrate
    /// \param[in] obj specific parameters object
    /// \param[in] range (optional) integration range
    template <typename T>
    inline double integrate(const utils::Function1D& func, const T& obj, const Limits& range = {}) const {
      return integrate(func, (void*)&obj, range);
    }
    /// Evaluate the integral of a function at a given value
    /// \param[in] func function to integrate
    /// \param[in] obj (optional) parameters object
    /// \param[in] range (optional) integration range
    virtual double integrate(const utils::Function1D& func, void* obj = nullptr, const Limits& range = {}) const = 0;

  protected:
    const Limits range_;
    const ParametersList func_params_;
    const int verbosity_;  ///< Integrator verbosity
  };
}  // namespace cepgen

#endif
