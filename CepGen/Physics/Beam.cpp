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
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  const double Beam::kMinKTFlux = 1.e-20;

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
      default:  // mode is invalid; try to infer it from KT flux
        switch (kt_flux_) {
          case KTFlux::P_Photon_Elastic:
          case KTFlux::P_Photon_Elastic_Budnev:
            mode_ = Mode::ProtonElastic;
            break;
          case KTFlux::P_Photon_Inelastic:
          case KTFlux::P_Photon_Inelastic_Budnev:
            mode_ = Mode::ProtonInelastic;
            break;
          case KTFlux::HI_Photon_Elastic:
          default:
            if (pdg_ == PDG::electron)
              mode_ = Mode::PointLikeFermion;
            else
              mode_ = Mode::Other;
            break;
        }
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

  pdgid_t Beam::daughterId() const {
    switch (kt_flux_) {
      case KTFlux::P_Gluon_KMR:
        return PDG::gluon;
      case KTFlux::P_Photon_Elastic:
      case KTFlux::P_Photon_Elastic_Budnev:
      case KTFlux::P_Photon_Inelastic:
      case KTFlux::P_Photon_Inelastic_Budnev:
      case KTFlux::HI_Photon_Elastic:
        return PDG::photon;
      case KTFlux::invalid:
      default:
        return PDG::invalid;
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
      flux = ktFluxHI(kt_flux_, x, q2, HeavyIon(pdg_));
    else {
      if (mx2 < 0.)
        throw CG_FATAL("Beam") << "Diffractive mass squared mX^2 should be specified!";
      if (!ff)
        throw CG_FATAL("Beam") << "Form factors should be specified!";
      if (!sf)
        throw CG_FATAL("Beam") << "Inelastic structure functions should be specified!";
      const auto mi = PDG::get().mass(pdg_), mi2 = mi * mi;
      flux = ktFluxNucl(kt_flux_, x, q2, ff, sf, mi2, mx2);
    }
    CG_DEBUG_LOOP("Beam:flux") << "Flux for (x=" << x << ", kT^2=" << q2 << "): " << flux << ".";
    return flux;
  }

  ParametersDescription Beam::description() {
    auto desc = ParametersDescription();
    desc.add<int>("pdgId", (int)PDG::proton);
    desc.add<double>("pz", 6500.);
    desc.add<int>("mode", (int)Beam::Mode::invalid);
    desc.add<int>("ktFlux", (int)KTFlux::invalid);
    return desc;
  }

  double Beam::ktFluxNucl(const KTFlux& type,
                          double x,
                          double kt2,
                          formfac::Parameterisation* ff,
                          strfun::Parameterisation* sf,
                          double mi2,
                          double mf2) {
    if (mi2 < 0.)
      mi2 = PDG::get().mass(PDG::proton);
    switch (type) {
      case KTFlux::P_Photon_Elastic:
      case KTFlux::P_Photon_Elastic_Budnev: {
        if (!ff)
          throw CG_FATAL("Beam:ktFlux") << "Elastic kT flux requires a modelling of electromagnetic form factors!";
        const double x2 = x * x;
        const double q2min = x2 * mi2 / (1. - x), q2 = q2min + kt2 / (1. - x);
        const double qnorm = 1. - q2min / q2;
        const auto& formfac = (*ff)(Beam::Mode::ProtonElastic, q2);
        const auto prefac = constants::ALPHA_EM * M_1_PI / q2;
        if (type == KTFlux::P_Photon_Elastic)
          return prefac * formfac.FE * qnorm * qnorm;
        const double f_D = formfac.FE * (1. - x) * qnorm;
        const double f_C = formfac.FM;
        return prefac * (1. - x) * (f_D + 0.5 * x2 * f_C);
      }
      case KTFlux::P_Photon_Inelastic:
      case KTFlux::P_Photon_Inelastic_Budnev: {
        if (!sf)
          throw CG_FATAL("Beam:ktFlux") << "Inelastic kT flux requires a modelling of structure functions!";
        const double x2 = x * x;
        const double q2min = (x * (mf2 - mi2) + x2 * mi2) / (1. - x), q2 = q2min + kt2 / (1. - x);
        const double qnorm = 1. - q2min / q2;
        const double denom = 1. / (q2 + mf2 - mi2), xbj = denom * q2;
        const auto prefac = constants::ALPHA_EM * M_1_PI * (1. - x) / q2;
        if (type == KTFlux::P_Photon_Inelastic)
          return prefac * sf->F2(xbj, q2) * denom * qnorm * qnorm;
        const double f_D = sf->F2(xbj, q2) * denom * (1. - x) * qnorm;
        const double f_C = sf->F1(xbj, q2) * 2. / q2;
        return prefac * (f_D + 0.5 * x2 * f_C);
      }
      case KTFlux::P_Gluon_KMR: {
        return kmr::GluonGrid::get()(x, kt2, mf2);
      } break;
      default:
        throw CG_FATAL("KTFlux") << "Invalid flux type: " << type;
    }
  }

  double Beam::ktFluxHI(const KTFlux& type, double x, double kt2, const HeavyIon& hi) {
    double flux = 0.;
    switch (type) {
      case KTFlux::HI_Photon_Elastic: {
        const auto mp = PDG::get().mass(PDG::proton);
        const double r_a = 1.1 * cbrt(hi.A), a0 = 0.7, m_a = hi.A * mp;
        const double q2_ela = (kt2 + x * x * m_a * m_a) / (1. - x),
                     cons = sqrt(q2_ela) / (constants::GEVM1_TO_M * 1e15);
        const double tau = cons * r_a, tau1 = cons * a0;
        // "Realistic nuclear form-factor" as used in STARLIGHT
        const double ff1 = 3. * (sin(tau) - tau * cos(tau)) / pow(tau + 1.e-10, 3);
        const double ff2 = 1. / (1. + tau1 * tau1);
        const double ela1 = pow(kt2 / (kt2 + x * x * m_a * m_a), 2);
        const double ela2 = pow(ff1 * ff2, 2) /*, ela3 = 1.-( q2_ela-kt2 )/q2_ela*/;
        const unsigned int z = (unsigned short)hi.Z;
        flux = constants::ALPHA_EM * M_1_PI * z * z * ela1 * ela2 / q2_ela;
      } break;
      default:
        throw CG_FATAL("KTFlux") << "Invalid flux type: " << type;
    }
    if (flux < kMinKTFlux)
      return 0.;
    return flux;
  }

  std::ostream& operator<<(std::ostream& os, const Beam& beam) {
    if (HeavyIon::isHI(beam.pdg_))
      os << HeavyIon(beam.pdg_);
    else
      os << (PDG::Id)beam.pdg_;
    os << " (" << beam.momentum_.pz() << " GeV/c), " << beam.mode_;
    if (beam.kt_flux_ != Beam::KTFlux::invalid)
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
      case Beam::Mode::Other:
        return os << "other";
    }
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const Beam::KTFlux& type) {
    switch (type) {
      case Beam::KTFlux::P_Photon_Elastic:
        return os << "elastic photon from proton";
      case Beam::KTFlux::P_Photon_Elastic_Budnev:
        return os << "elastic photon from proton (Budnev)";
      case Beam::KTFlux::P_Photon_Inelastic:
        return os << "inelastic photon from proton";
      case Beam::KTFlux::P_Photon_Inelastic_Budnev:
        return os << "inelastic photon from proton (Budnev)";
      case Beam::KTFlux::P_Gluon_KMR:
        return os << "elastic gluon from proton (KMR)";
      case Beam::KTFlux::HI_Photon_Elastic:
        return os << "elastic photon from HI";
      case Beam::KTFlux::invalid:
      default:
        return os << "unrecognised flux (" << (int)type << ")";
    }
  }
}  // namespace cepgen
