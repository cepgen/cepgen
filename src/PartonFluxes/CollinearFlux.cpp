/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

#include "CepGen/PartonFluxes/CollinearFlux.h"
#include "CepGen/Physics/Utils.h"

namespace cepgen {
  CollinearFlux::CollinearFlux(const ParametersList& params) : PartonFlux(params) {}

  ParametersDescription CollinearFlux::description() {
    auto desc = PartonFlux::description();
    desc.setDescription("Collinear parton flux");
    return desc;
  }

  double CollinearFlux::fluxQ2(double x, double q2) const {
    return fluxMX2(x, utils::mX2(x, q2, mass2()) /*FIXME on-shell assumption: xbj == x*/);
  }

  double CollinearFlux::fluxMX2(double x, double mx2) const {
    return fluxQ2(x, utils::q2(x, mass2(), mx2 <= 0. ? mass2() : mx2) /*FIXME on-shell assumption: xbj == x*/);
  }
}  // namespace cepgen
