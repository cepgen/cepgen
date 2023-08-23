/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
        pdg_id_(steer<pdgid_t>("pdgId")),
        momentum_(Momentum::fromPxPyPzM(
            0.,
            0.,
            steer<double>("pz"),
            HeavyIon::isHI(pdg_id_) ? HeavyIon::mass(HeavyIon::fromPdgId(pdg_id_)) : PDG::get().mass(pdg_id_))),
        mode_(steerAs<int, Mode>("mode")) {
    if (pdg_id_ == PDG::electron)
      mode_ = Mode::PointLikeFermion;
    else if (HeavyIon::isHI(pdg_id_))
      mode_ = Mode::HIElastic;
  }

  void Beam::initialise() {
    const auto flux_info = steer<ParametersList>("partonFlux").set("pdgId", pdg_id_);
    if (!flux_info.name<std::string>().empty())
      flux_ = PartonFluxFactory::get().build(flux_info);
  }

  bool Beam::fragmented() const {
    if (flux_)
      return flux().fragmenting();
    return mode_ != Mode::HIElastic && mode_ != Mode::ProtonElastic;
  }

  pdgid_t Beam::daughterId() const { return flux().partonPdgId(); }

  const PartonFlux& Beam::flux() const {
    if (!flux_)
      throw CG_FATAL("Beam:flux") << "Beam flux requested although it was not yet initialised.";
    return *flux_;
  }

  ParametersDescription Beam::description() {
    auto desc = ParametersDescription();
    desc.addAs<int, pdgid_t>("pdgId", PDG::proton);
    desc.add<double>("pz", 0.);
    desc.addAs<int, Beam::Mode>("mode", Beam::Mode::invalid);
    desc.add<ParametersDescription>("partonFlux", ParametersDescription());
    return desc;
  }

  std::ostream& operator<<(std::ostream& os, const Beam& beam) {
    if (HeavyIon::isHI(beam.pdg_id_))
      os << HeavyIon::fromPdgId(beam.pdg_id_);
    else
      os << (PDG::Id)beam.pdg_id_;
    os << " (" << beam.momentum_.pz() << " GeV/c), " << beam.mode_;
    if (beam.flux_)
      os << " [part.flux: " << beam.flux_->name() << "]";
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const Beam::Mode& mode) {
    switch (mode) {
      case Beam::Mode::invalid:
        return os << "{invalid}";
      case Beam::Mode::ProtonElastic:
        return os << "el.proton";
      case Beam::Mode::HIElastic:
        return os << "el.ion";
      case Beam::Mode::PointLikeScalar:
        return os << "gen.scalar";
      case Beam::Mode::PointLikeFermion:
        return os << "gen.fermion";
      case Beam::Mode::CompositeScalar:
        return os << "comp.scalar";
      case Beam::Mode::ProtonInelastic:
        return os << "inel.proton";
      case Beam::Mode::Other:
        return os << "other";
    }
    return os;
  }
}  // namespace cepgen
