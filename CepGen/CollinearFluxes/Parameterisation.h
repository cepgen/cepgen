/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2023  Laurent Forthomme
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

#include "CepGen/PartonFluxes/PartonFlux.h"

namespace cepgen {
  /// Collinear fluxes modelling scope
  namespace collflux {
    /// Generic collinear flux parameterisation
    class Parameterisation : public PartonFlux {
    public:
      /// User-steered parameterisation object constructor
      explicit Parameterisation(const ParametersList&);

      /// Generic description for the collinear flux
      static ParametersDescription description();

      virtual double operator()(double, double = 0.) const = 0;
      double operator()(double /*x*/, double /*kt2*/, double /*mf2*/) const override final;
      bool ktFactorised() const final { return false; }

    protected:
      Limits q2_range_;
      const double qscale_{0.};
    };
  }  // namespace collflux
}  // namespace cepgen

#endif
