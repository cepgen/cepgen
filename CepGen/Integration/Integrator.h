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

#ifndef CepGen_Integration_Integrator_h
#define CepGen_Integration_Integrator_h

#include <functional>

#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Utils/Limits.h"
#include "CepGen/Utils/Value.h"

namespace cepgen {
  class Integrand;
}  // namespace cepgen

namespace cepgen {
  /// Integration algorithm
  class Integrator : public NamedModule<Integrator> {
  public:
    explicit Integrator(const ParametersList&);

    static ParametersDescription description();

    virtual double eval(Integrand&, const std::vector<double>&) const;  ///< Compute function value at one point

    /// Evaluate the integral for a given range
    Value integrate(Integrand& integrand, const std::vector<Limits>& = {});
    /// Evaluate the integral of a function for a given range
    /// \param[in] integrand Function to integrate
    /// \param[in] range_1d integration range
    Value integrate(const std::function<double(double)>& integrand, const Limits& range_1d = {0., 1.});
    /// Evaluate the integral of a function for a given range
    /// \param[in] integrand Function to integrate
    /// \param[in] range integration range
    Value integrate(const std::function<double(const std::vector<double>&)>& integrand,
                    const std::vector<Limits>& range);

  protected:
    /// Evaluate the integral of a function for a given range
    /// \param[in] integrand Function to integrate
    /// \param[in] range integration range
    virtual Value run(Integrand& integrand, const std::vector<Limits>& range) = 0;

    const int verbosity_;  ///< Integrator verbosity
  };
}  // namespace cepgen

#endif
