/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  Beam::Beam(const ParametersList& params)
      : SteeredObject(params),
        pdg_(steerAs<int, pdgid_t>("pdgId")),
        momentum_(Momentum::fromPxPyPzM(
            0., 0., steer<double>("pz"), HeavyIon::isHI(pdg_) ? HeavyIon::mass(HeavyIon(pdg_)) : PDG::get().mass(pdg_))),
        mode_(steerAs<int, Mode>("mode")),
        kt_flux_(steerAs<int, KTFlux>("ktFlux")) {
    if (pdg_ == PDG::electron)
      mode_ = Mode::PointLikeFermion;
    switch (mode_) {
      case Mode::ProtonElastic:
        if (kt_flux_ != KTFlux::P_Photon_Elastic && kt_flux_ != KTFlux::P_Photon_Elastic_Budnev &&
            kt_flux_ != KTFlux::HI_Photon_Elastic && kt_flux_ != KTFlux::P_Gluon_KMR) {
          kt_flux_ = HeavyIon::isHI(pdg_) ? KTFlux::HI_Photon_Elastic : KTFlux::P_Photon_Elastic_Budnev;
          CG_DEBUG("Beam") << "KT flux for incoming parton set to \"" << kt_flux_ << "\".";
        }
        break;
      case Mode::ProtonInelastic:
        if (kt_flux_ != KTFlux::P_Photon_Inelastic && kt_flux_ != KTFlux::P_Photon_Inelastic_Budnev) {
          if (HeavyIon::isHI(pdg_))
            throw CG_FATAL("KTProcess:kinematics") << "Inelastic photon emission from HI not yet supported!";
          kt_flux_ = KTFlux::P_Photon_Inelastic_Budnev;
          CG_DEBUG("Beam") << "KT flux for incoming parton set to \"" << kt_flux_ << "\".";
        }
        break;
      default:
        break;
    }
  }

  bool Beam::fragmented() const {
    switch (mode_) {
      case Mode::ProtonElastic:
      case Mode::PointLikeScalar:
      case Mode::PointLikeFermion:
      default:
        return false;
      case Mode::ProtonInelastic:
        return true;
    }
  }

  Beam::FormFactors Beam::flux(double q2,
                               double mx2,
                               formfac::Parameterisation* ff,
                               strfun::Parameterisation* sf) const {
    const auto flux = (*ff)(mode_, q2, mx2, sf);
    CG_DEBUG_LOOP("Beam:flux") << "Flux for (q2=" << q2 << ", mX2=" << mx2 << "): "
                               << "FE=" << flux.FE << ", FM=" << flux.FM << ".";
    return FormFactors{flux.FE, flux.FM};
  }

  double Beam::ktFlux(
      double x, double q2, double mx2, formfac::Parameterisation* ff, strfun::Parameterisation* sf) const {
    double flux = 0.;
    if (HeavyIon::isHI(pdg_))  // check if we are in heavy ion mode
      flux = cepgen::ktFlux((KTFlux)kt_flux_, x, q2, HeavyIon(pdg_));
    else {
      if (mx2 < 0.)
        throw CG_FATAL("Beam") << "Diffractive mass squared mX^2 should be specified!";
      if (!ff)
        throw CG_FATAL("Beam") << "Form factors should be specified!";
      if (!sf)
        throw CG_FATAL("Beam") << "Inelastic structure functions should be specified!";
      const auto mi = PDG::get().mass(pdg_), mi2 = mi * mi;
      flux = cepgen::ktFlux((KTFlux)kt_flux_, x, q2, *ff, *sf, mi2, mx2);
    }
    CG_DEBUG_LOOP("Beam:flux") << "Flux for (x=" << x << ", kT^2=" << q2 << "): " << flux << ".";
    return flux;
  }

  ParametersDescription Beam::description() {
    auto desc = ParametersDescription();
    desc.add<int>("pdgId", (int)PDG::proton);
    desc.add<int>("pz", 6500.);
    desc.add<int>("mode", (int)Beam::Mode::invalid);
    desc.add<int>("ktFlux", (int)KTFlux::invalid);
    return desc;
  }

  std::ostream& operator<<(std::ostream& os, const Beam& beam) {
    if (HeavyIon::isHI(beam.pdg_))
      os << HeavyIon(beam.pdg_);
    else
      os << (PDG::Id)beam.pdg_;
    os << " (" << beam.momentum_.pz() << " GeV/c), " << beam.mode_;
    if (beam.kt_flux_ != KTFlux::invalid)
      os << " [unint.flux: " << beam.kt_flux_ << "]";
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const Beam::Mode& mode) {
    switch (mode) {
      case Beam::Mode::invalid:
        return os << "{invalid}";
      case Beam::Mode::ProtonElastic:
        return os << "el.proton";
      case Beam::Mode::PointLikeScalar:
        return os << "gen.scalar";
      case Beam::Mode::PointLikeFermion:
        return os << "gen.fermion";
      case Beam::Mode::CompositeScalar:
        return os << "comp.scalar";
      case Beam::Mode::ProtonInelastic:
        return os << "inel.proton";
    }
    return os;
  }
}  // namespace cepgen
