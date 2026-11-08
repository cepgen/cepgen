/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#ifndef CepGen_PartonFluxes_CollinearFlux_h
#define CepGen_PartonFluxes_CollinearFlux_h

#include "CepGen/PartonFluxes/PartonFlux.h"

namespace cepgen {
  /// Base object for a collinear parton flux parameterisation
  class CollinearFlux : public PartonFlux {
  public:
    explicit CollinearFlux(const ParametersList&);

    static ParametersDescription description();

    virtual double fluxQ2(double x, double q2) const;         ///< Compute the collinear flux for an x/virtuality
    virtual double fluxMX2(double x, double mf2 = 0.) const;  ///< Compute the collinear flux for an x/remnant mass

    bool ktFactorised() const final { return false; }
  };
}  // namespace cepgen

#endif
