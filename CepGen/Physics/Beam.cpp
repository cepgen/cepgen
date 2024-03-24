/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  Beam::Beam(const ParametersList& params)
      : SteeredObject(params),
        pdg_id_(steerAs<int, pdgid_t>("pdgId")),
        momentum_(Momentum::fromPxPyPzM(0., 0., steer<double>("pz"), PDG::get().mass(pdg_id_))),
        flux_info_(steer<ParametersList>("partonFlux")),
        elastic_(steer<bool>("elastic")) {}

  ParametersDescription Beam::description() {
    auto desc = ParametersDescription();
    desc.addAs<int, pdgid_t>("pdgId", PDG::proton);
    desc.add<double>("pz", 0.);
    desc.add<ParametersDescription>("partonFlux", ParametersDescription());
    return desc;
  }

  std::ostream& operator<<(std::ostream& os, const Beam& beam) {
    if (HeavyIon::isHI(beam.pdg_id_))
      os << HeavyIon::fromPdgId(beam.pdg_id_);
    else
      os << (PDG::Id)beam.pdg_id_;
    os << " (" << beam.momentum_.pz() << " GeV/c) " << (beam.elastic_ ? "elastic" : "inelastic");
    if (!beam.flux_info_.name<std::string>().empty())
      os << " [part.flux: " << beam.flux_info_ << "]";
    return os;
  }
}  // namespace cepgen
