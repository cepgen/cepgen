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

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/PartonFlux.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/KTPhaseSpaceGenerator.h"

namespace cepgen {
  namespace proc {
    FactorisedProcess::FactorisedProcess(const ParametersList& params, const pdgids_t& central)
        : Process(params), produced_parts_(central), psgen_(new KTPhaseSpaceGenerator(this)) {
      event().map()[Particle::CentralSystem].resize(central.size());
    }

    FactorisedProcess::FactorisedProcess(const FactorisedProcess& proc)
        : Process(proc), produced_parts_(proc.produced_parts_), psgen_(new KTPhaseSpaceGenerator(this)) {}

    void FactorisedProcess::addEventContent() {
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
      CG_DEBUG("FactorisedProcess:addEventContent")
          << "Addition of:\n\t"
          << "Intermediate partons: "
          << pdgids_t{psgen_->positiveFlux().partonPdgId(), psgen_->negativeFlux().partonPdgId()} << "\n\t"
          << "Produced system: " << produced_parts_ << ".\n\t" << event();
    }

    void FactorisedProcess::prepareKinematics() {
      if (!psgen_)
        throw CG_FATAL("FactorisedProcess:prepareKinematics")
            << "Phase space generator not set. Please check your process initialisation procedure, as you might "
               "be doing something irregular.";
      psgen_->init();

      event().oneWithRole(Particle::Parton1).setPdgId(psgen_->positiveFlux().partonPdgId());
      event().oneWithRole(Particle::Parton2).setPdgId(psgen_->negativeFlux().partonPdgId());

      //============================================================================================
      // register all process-dependent variables
      //============================================================================================

      preparePhaseSpace();

      //============================================================================================
      // register the outgoing remnants' variables
      //============================================================================================

      mX2() = pA().mass2();
      if (!kinematics().incomingBeams().positive().elastic())
        defineVariable(
            mX2(), Mapping::square, kinematics().cuts().remnants.mx, "Positive z proton remnant squared mass");

      mY2() = pB().mass2();
      if (!kinematics().incomingBeams().negative().elastic())
        defineVariable(
            mY2(), Mapping::square, kinematics().cuts().remnants.mx, "Negative z proton remnant squared mass");
    }

    double FactorisedProcess::computeWeight() {
      if (const auto cent_me = computeFactorisedMatrixElement() > 0.)
        return psgen_->fluxes() * cent_me;
      return 0.;  // avoid computing the fluxes if the matrix element is already null
    }

    void FactorisedProcess::fillKinematics(bool) {
      fillCentralParticlesKinematics();  // process-dependent!

      // beam systems
      if (!kinematics().incomingBeams().positive().elastic())
        pX().setMass2(mX2());
      if (!kinematics().incomingBeams().negative().elastic())
        pY().setMass2(mY2());

      // parton systems
      auto& part1 = event().oneWithRole(Particle::Parton1);
      auto& part2 = event().oneWithRole(Particle::Parton2);
      part1.setMomentum(event().oneWithRole(Particle::IncomingBeam1).momentum() - pX(), true);
      part2.setMomentum(event().oneWithRole(Particle::IncomingBeam2).momentum() - pY(), true);

      // two-parton system
      event().oneWithRole(Particle::Intermediate).setMomentum(part1.momentum() + part2.momentum());
    }

    ParametersDescription FactorisedProcess::description() {
      auto desc = Process::description();
      desc.setDescription("Unnamed factorised process");
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
