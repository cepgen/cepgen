/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#ifndef CepGen_CollinearFluxes_CollinearFlux_h
#define CepGen_CollinearFluxes_CollinearFlux_h

#include "CepGen/Physics/PartonFlux.h"

namespace cepgen {
  class CollinearFlux : public PartonFlux {
  public:
    explicit CollinearFlux(const ParametersList&);

    static ParametersDescription description();

    /// Compute the collinear flux for this x value and virtuality
    virtual double fluxQ2(double x, double q2) const;
    /// Compute the collinear flux for this x value and remnant mass
    virtual double fluxMX2(double x, double mf2 = 0.) const;

    bool ktFactorised() const override final { return false; }
  };
}  // namespace cepgen

#endif
