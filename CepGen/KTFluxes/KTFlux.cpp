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

#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Physics/PartonFlux.h"

namespace cepgen {
  KTFlux::KTFlux(const ParametersList& params) : PartonFlux(params) {}

  ParametersDescription KTFlux::description() {
    auto desc = PartonFlux::description();
    desc.setDescription("kT-factorised flux");
    return desc;
  }

  double KTFlux::fluxQ2(double x, double kt2, double q2) const { return fluxMX2(x, kt2, mX2(x, kt2, q2)); }

  double KTFlux::fluxMX2(double x, double kt2, double mf2) const { return fluxQ2(x, kt2, computeQ2(x, kt2, mf2).q2); }

  KTFlux::Q2Values KTFlux::computeQ2(double x, double kt2, double mx2) const {
    Q2Values out;
    const auto mi2 = mass2();
    const auto dm2 = (mx2 == 0.) ? 0. : mx2 - mi2;
    out.min = ((x * dm2) + x * x * mi2) / (1. - x);
    out.q2 = out.min + kt2 / (1. - x);
    return out;
  }

  double KTFlux::mX2(double x, double kt2, double q2) const {
    return mass2() + (q2 * (1. - x) - kt2 - x * x * mass2()) / x;
  }
}  // namespace cepgen
