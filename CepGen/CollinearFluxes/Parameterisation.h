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

#ifndef CepGen_CollinearFluxes_Parameterisation_h
#define CepGen_CollinearFluxes_Parameterisation_h

#include <iosfwd>
#include <memory>

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  /// Collinear fluxes modelling scope
  namespace collflux {
    /// Generic collinear flux parameterisation
    class Parameterisation : public NamedModule<std::string> {
    public:
      /// User-steered parameterisation object constructor
      explicit Parameterisation(const ParametersList&);
      virtual ~Parameterisation() = default;

      /// Generic description for the structure functions
      static ParametersDescription description();

      /// Compute the collinear flux for this x value
      virtual double operator()(double /*x*/, double /*mx*/ = 0.) const { return 0.; }

    protected:
      const double mp_{0.};   ///< Proton mass, in GeV/c^2
      const double mp2_{0.};  ///< Squared proton mass, in GeV^2/c^4
      Limits q2_range_;
      const double qscale_{0.};
    };
  }  // namespace collflux
}  // namespace cepgen

#endif
