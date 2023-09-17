/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2023  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/Process/KTProcess.h"

namespace cepgen {
  namespace proc {
    KTProcess::KTProcess(const ParametersList& params, const pdgids_t& central)
        : Process(params), produced_parts_(central) {
      event().map()[Particle::CentralSystem].resize(central.size());
    }

    void KTProcess::addEventContent() {
      Process::setEventContent(
          {// incoming state
           {Particle::IncomingBeam1, kinematics().incomingBeams().positive().pdgId()},
           {Particle::IncomingBeam2, kinematics().incomingBeams().negative().pdgId()},
           {Particle::Parton1, PDG::invalid},
           {Particle::Parton2, PDG::invalid}},
          {// outgoing state
           {Particle::OutgoingBeam1, {kinematics().incomingBeams().positive().pdgId()}},
           {Particle::OutgoingBeam2, {kinematics().incomingBeams().negative().pdgId()}},
           {Particle::CentralSystem, produced_parts_}});
      setExtraContent();
    }

    void KTProcess::prepareKinematics() {
      auto set_beam_properties = [](Beam& beam, std::shared_ptr<PartonFlux>& flux) {
        auto params = beam.partonFluxParameters();
        if (beam.elastic()) {
          if (HeavyIon::isHI(beam.pdgId()))
            params = PartonFluxFactory::get().describeParameters("ElasticHeavyIonKT").validate(params);
          else
            params = PartonFluxFactory::get().describeParameters("BudnevElasticKT").validate(params);
        } else
          params = PartonFluxFactory::get().describeParameters("BudnevInelasticKT").validate(params);
        flux.reset(PartonFluxFactory::get().build(params).release());
      };
      set_beam_properties(kinematics().incomingBeams().positive(), pos_flux_);
      set_beam_properties(kinematics().incomingBeams().negative(), neg_flux_);

      if (!pos_flux_ || !neg_flux_)
        throw CG_FATAL("KTProcess:prepareKinematics") << "Invalid incoming parton fluxes.";

      event().oneWithRole(Particle::Parton1).setPdgId(pos_flux_->partonPdgId());
      event().oneWithRole(Particle::Parton2).setPdgId(neg_flux_->partonPdgId());

      CG_DEBUG("KTProcess:prepareKinematics")
          << "Partons: " << pdgids_t{pos_flux_->partonPdgId(), neg_flux_->partonPdgId()} << ", "
          << "central system: " << produced_parts_ << ". " << event();

      //============================================================================================
      // register the incoming partons' variables
      //============================================================================================

      const auto log_lim_kt = kinematics().cuts().initial.qt.compute(std::log).truncate(Limits{-10., 10.});
      defineVariable(m_qt1_, Mapping::exponential, log_lim_kt, "Positive-z parton virtuality");
      defineVariable(m_qt2_, Mapping::exponential, log_lim_kt, "Negative-z parton virtuality");

      const auto lim_phi = kinematics().cuts().initial.phi_qt.truncate(Limits{0., 2. * M_PI});
      defineVariable(m_phi_qt1_, Mapping::linear, lim_phi, "Positive-z parton azimuthal angle");
      defineVariable(m_phi_qt2_, Mapping::linear, lim_phi, "Negative-z parton azimuthal angle");

      //============================================================================================
      // register all process-dependent variables
      //============================================================================================

      preparePhaseSpace();

      //============================================================================================
      // register the outgoing remnants' variables
      //============================================================================================

      mX2() = pA().mass2();
      if (!kinematics().incomingBeams().positive().elastic())
        defineVariable(mX2(), Mapping::square, kinematics().cuts().remnants.mx, "Positive z-beam remnant squared mass");
      mY2() = pB().mass2();
      if (!kinematics().incomingBeams().negative().elastic())
        defineVariable(mY2(), Mapping::square, kinematics().cuts().remnants.mx, "Negative z-beam remnant squared mass");
    }

    double KTProcess::computeWeight() {
      // compute the transverse kinematics of the initial partons
      q1() = Momentum::fromPtEtaPhiE(m_qt1_, 0., m_phi_qt1_);
      q2() = Momentum::fromPtEtaPhiE(m_qt2_, 0., m_phi_qt2_);

      // compute the central matrix element
      const auto cent_me = computeKTFactorisedMatrixElement();
      if (cent_me <= 0)
        return 0.;  // avoid computing the fluxes if the matrix element is already null

      // convolute with fluxes according to modelling specified in parameters card
      const auto& flux1 = dynamic_cast<const KTFlux&>(*pos_flux_);
      const auto& flux2 = dynamic_cast<const KTFlux&>(*neg_flux_);

      return (flux1.fluxMX2(x1_, m_qt1_ * m_qt1_, mX2()) * M_1_PI) *
             (flux2.fluxMX2(x2_, m_qt2_ * m_qt2_, mY2()) * M_1_PI) * cent_me;
    }

    void KTProcess::fillKinematics(bool) {
      t1() = utils::kt::q2(x1_, m_qt1_ * m_qt1_, mA2(), mX2());
      t2() = utils::kt::q2(x2_, m_qt2_ * m_qt2_, mB2(), mY2());

      fillCentralParticlesKinematics();  // process-dependent!

      // beam systems
      if (!kinematics().incomingBeams().positive().elastic())
        pX().setMass2(mX2());
      if (!kinematics().incomingBeams().negative().elastic())
        pY().setMass2(mY2());

      // parton systems
      auto& p1 = event().oneWithRole(Particle::Parton1);
      auto& p2 = event().oneWithRole(Particle::Parton2);
      p1.setMomentum(pA() - pX(), true);
      p2.setMomentum(pB() - pY(), true);

      // two-parton system
      event().oneWithRole(Particle::Intermediate).setMomentum(p1.momentum() + p2.momentum(), true);
    }

    ParametersDescription KTProcess::description() {
      auto desc = Process::description();
      desc.setDescription("Unnamed kT-factorised process");
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
