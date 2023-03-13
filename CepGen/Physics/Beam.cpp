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
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  Beam::Beam(const ParametersList& params)
      : SteeredObject(params),
        pdg_(steerAs<int, pdgid_t>("pdgId")),
        momentum_(Momentum::fromPxPyPzM(
            0.,
            0.,
            steer<double>("pz"),
            HeavyIon::isHI(pdg_) ? HeavyIon::mass(HeavyIon::fromPdgId(pdg_)) : PDG::get().mass(pdg_))),
        mode_(steerAs<int, Mode>("mode")) {
    if (pdg_ == PDG::electron)
      mode_ = Mode::PointLikeFermion;
    else if (HeavyIon::isHI(pdg_))
      mode_ = Mode::HIElastic;
  }

  void Beam::initialise() {
    const auto& pflux_params = steer<ParametersList>("partonFlux");
    switch (mode_) {
      case Mode::HIElastic:
        flux_ = PartonFluxFactory::get().build("ElasticHeavyIonKT", params_ + pflux_params);
        break;
      case Mode::ProtonInelastic:
        flux_ = PartonFluxFactory::get().build("BudnevInelasticKT", params_ + pflux_params);
        break;
      case Mode::ProtonElastic:
        flux_ = PartonFluxFactory::get().build("BudnevElasticKT", params_ + pflux_params);
        break;
      default:
        CG_WARNING("Beam:initialise") << "Invalid beam mode retrieved: '" << mode_ << "'.";
        break;
    }
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

  double Beam::flux(double x, double q2, double mx2) const { return flux()(x, q2, mx2); }

  ParametersDescription Beam::description() {
    auto desc = ParametersDescription();
    desc.addAs<int, pdgid_t>("pdgId", PDG::proton);
    desc.add<double>("pz", 0.);
    desc.add<int>("mode", (int)Beam::Mode::invalid);
    desc.add<ParametersDescription>("partonFlux", PartonFluxFactory::get().describeParameters("BudnevElasticKT"));
    return desc;
  }

  std::ostream& operator<<(std::ostream& os, const Beam& beam) {
    if (HeavyIon::isHI(beam.pdg_))
      os << HeavyIon::fromPdgId(beam.pdg_);
    else
      os << (PDG::Id)beam.pdg_;
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
