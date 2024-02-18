/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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
#include "CepGen/Process/CollinearPhaseSpaceGenerator.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/KTPhaseSpaceGenerator.h"
#include "CepGen/Utils/Math.h"

namespace cepgen {
  namespace proc {
    FactorisedProcess::FactorisedProcess(const ParametersList& params, const pdgids_t& central)
        : Process(params),
          produced_parts_(central),
          psgen_(steer<bool>("ktFactorised")
                     ? std::unique_ptr<PhaseSpaceGenerator>(new KTPhaseSpaceGenerator(this))
                     : std::unique_ptr<PhaseSpaceGenerator>(new CollinearPhaseSpaceGenerator(this))),
          store_alphas_(steer<bool>("storeAlphas")) {
      event().map()[Particle::CentralSystem].resize(central.size());
    }

    FactorisedProcess::FactorisedProcess(const FactorisedProcess& proc)
        : Process(proc),
          produced_parts_(proc.produced_parts_),
          psgen_(proc.psgen_->ktFactorised()
                     ? std::unique_ptr<PhaseSpaceGenerator>(new KTPhaseSpaceGenerator(this))
                     : std::unique_ptr<PhaseSpaceGenerator>(new CollinearPhaseSpaceGenerator(this))),
          store_alphas_(proc.store_alphas_) {}

    void FactorisedProcess::addEventContent() {
      Process::setEventContent({{Particle::IncomingBeam1, {kinematics().incomingBeams().positive().pdgId()}},
                                {Particle::IncomingBeam2, {kinematics().incomingBeams().negative().pdgId()}},
                                {Particle::OutgoingBeam1, {kinematics().incomingBeams().positive().pdgId()}},
                                {Particle::OutgoingBeam2, {kinematics().incomingBeams().negative().pdgId()}},
                                {Particle::CentralSystem, produced_parts_}});
    }

    void FactorisedProcess::prepareKinematics() {
      if (!psgen_)
        throw CG_FATAL("FactorisedProcess:prepareKinematics")
            << "Phase space generator not set. Please check your process initialisation procedure, as you might "
               "be doing something irregular.";
      psgen_->initialise();

      event().oneWithRole(Particle::Parton1).setPdgId(psgen_->positiveFlux().partonPdgId());
      event().oneWithRole(Particle::Parton2).setPdgId(psgen_->negativeFlux().partonPdgId());

      CG_DEBUG("FactorisedProcess:prepareKinematics")
          << "Partons: " << pdgids_t{psgen_->positiveFlux().partonPdgId(), psgen_->negativeFlux().partonPdgId()} << ", "
          << "central system: " << produced_parts_ << ". " << event();

      // register all process-dependent variables
      prepareFactorisedPhaseSpace();

      // register the outgoing remnants' variables
      if (!kinematics().incomingBeams().positive().elastic())
        defineVariable(mX2(), Mapping::square, kinematics().cuts().remnants.mx, "Positive-z beam remnant squared mass");
      if (!kinematics().incomingBeams().negative().elastic())
        defineVariable(mY2(), Mapping::square, kinematics().cuts().remnants.mx, "Negative-z beam remnant squared mass");
    }

    double FactorisedProcess::computeWeight() {
      if (!psgen_->generatePartonKinematics())
        return 0.;
      const auto cent_me = computeFactorisedMatrixElement();
      if (!utils::positive(cent_me))
        return 0.;  // avoid computing the fluxes if the matrix element is already null or invalid
      const auto fluxes_weight = psgen_->fluxes();
      if (!utils::positive(fluxes_weight))
        return 0.;
      return fluxes_weight * cent_me;
    }

    void FactorisedProcess::fillKinematics() {
      fillCentralParticlesKinematics();  // process-dependent!

      // beam systems
      if (!kinematics().incomingBeams().positive().elastic())
        pX().setMass2(mX2());
      if (!kinematics().incomingBeams().negative().elastic())
        pY().setMass2(mY2());

      // parton systems
      auto& part1 = event().oneWithRole(Particle::Parton1);
      auto& part2 = event().oneWithRole(Particle::Parton2);
      part1.setMomentum(pA() - pX(), true);
      part2.setMomentum(pB() - pY(), true);

      // add couplings to metadata
      if (store_alphas_) {
        const auto two_part_mass = (part1.momentum() + part2.momentum()).mass();
        event().metadata["alphaEM"] = alphaEM(two_part_mass);
        event().metadata["alphaS"] = alphaS(two_part_mass);
      }
    }

    ParametersDescription FactorisedProcess::description() {
      auto desc = Process::description();
      desc.setDescription("Unnamed factorised process");
      desc.add<bool>("storeAlphas", false)
          .setDescription("store the electromagnetic and strong coupling constants to the event content?");
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
