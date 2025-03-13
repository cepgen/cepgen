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

#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/PDG.h"

using namespace cepgen;

struct KMRGluonKTFlux final : KTFlux {
  using KTFlux::KTFlux;
  static ParametersDescription description() {
    auto desc = KTFlux::description();
    desc.setDescription("Proton inelastic gluon emission (KMR flux)");
    return desc;
  }
  double fluxMX2(double x, double kt2, double mx2) const override {
    if (!x_range_.contains(x, true))
      return 0.;
    return kmr::GluonGrid::get()(x, kt2, mx2);
  }
  pdgid_t partonPdgId() const override { return PDG::gluon; }
  bool fragmenting() const override { return false; }
  double mass2() const override { return mp2_; }
};
REGISTER_KT_FLUX("KMR", 20, KMRGluonKTFlux);
