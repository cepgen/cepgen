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

#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  struct KMRGluonKTFlux final : public KTFlux {
    using KTFlux::KTFlux;
    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Proton inelastic gluon emission (KMR flux)");
      return desc;
    }
    double operator()(double x, double kt2, double mx2) const override final {
      if (!x_range_.contains(x))
        return 0.;
      return kmr::GluonGrid::get()(x, kt2, mx2);
    }
    pdgid_t partonPdgId() const override final { return PDG::gluon; }
    bool fragmenting() const override { return false; }
  };
}  // namespace cepgen

REGISTER_FLUX("KMR", KMRGluonKTFlux);
