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
#include "CepGen/Modules/CentralPhaseSpaceGeneratorFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/PartonFlux.h"
#include "CepGen/Process/CentralPhaseSpaceGenerator.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsCollinearPhaseSpaceGenerator.h"
#include "CepGen/Process/PartonsKTPhaseSpaceGenerator.h"
#include "CepGen/Utils/Math.h"

namespace cepgen {
  namespace proc {
    FactorisedProcess::FactorisedProcess(const ParametersList& params, const pdgids_t& central)
        : Process(params),
          part_psgen_(steer<bool>("ktFactorised")
                          ? std::unique_ptr<PartonsPhaseSpaceGenerator>(new PartonsKTPhaseSpaceGenerator(this))
                          : std::unique_ptr<PartonsPhaseSpaceGenerator>(new PartonsCollinearPhaseSpaceGenerator(this))),
          cent_psgen_(CentralPhaseSpaceGeneratorFactory::get().build(steer<ParametersList>("kinematicsGenerator"))),
          store_alphas_(steer<bool>("storeAlphas")) {
      event().map()[Particle::CentralSystem].resize(central.size());
    }

    FactorisedProcess::FactorisedProcess(const FactorisedProcess& proc)
        : Process(proc),
          part_psgen_(proc.part_psgen_->ktFactorised()
                          ? std::unique_ptr<PartonsPhaseSpaceGenerator>(new PartonsKTPhaseSpaceGenerator(this))
                          : std::unique_ptr<PartonsPhaseSpaceGenerator>(new PartonsCollinearPhaseSpaceGenerator(this))),
          cent_psgen_(CentralPhaseSpaceGeneratorFactory::get().build(proc.cent_psgen_->parameters())),
          store_alphas_(proc.store_alphas_) {}

    void FactorisedProcess::addEventContent() {
      Process::setEventContent({{Particle::IncomingBeam1, {kinematics().incomingBeams().positive().pdgId()}},
                                {Particle::IncomingBeam2, {kinematics().incomingBeams().negative().pdgId()}},
                                {Particle::OutgoingBeam1, {kinematics().incomingBeams().positive().pdgId()}},
                                {Particle::OutgoingBeam2, {kinematics().incomingBeams().negative().pdgId()}},
                                {Particle::CentralSystem, cent_psgen_->particles()}});
    }

    void FactorisedProcess::prepareKinematics() {
      if (!part_psgen_)
        throw CG_FATAL("FactorisedProcess:prepareKinematics")
            << "Phase space generator not set. Please check your process initialisation procedure, as you might "
               "be doing something irregular.";
      part_psgen_->initialise();
      cent_psgen_->initialise(this);

      event().oneWithRole(Particle::Parton1).setPdgId(part_psgen_->positiveFlux().partonPdgId());
      event().oneWithRole(Particle::Parton2).setPdgId(part_psgen_->negativeFlux().partonPdgId());

      CG_DEBUG("FactorisedProcess:prepareKinematics")
          << "Partons: "
          << pdgids_t{part_psgen_->positiveFlux().partonPdgId(), part_psgen_->negativeFlux().partonPdgId()} << ", "
          << "central system: " << cent_psgen_->particles() << ". " << event();

      // register all process-dependent variables
      prepareFactorisedPhaseSpace();

      // register the outgoing remnants' variables
      if (!kinematics().incomingBeams().positive().elastic())
        defineVariable(mX2(), Mapping::square, kinematics().cuts().remnants.mx, "Positive-z beam remnant squared mass");
      if (!kinematics().incomingBeams().negative().elastic())
        defineVariable(mY2(), Mapping::square, kinematics().cuts().remnants.mx, "Negative-z beam remnant squared mass");
    }

    double FactorisedProcess::computeWeight() {
      if (!part_psgen_->generatePartonKinematics())
        return 0.;
      const auto cent_kin_weight = cent_psgen_->generateKinematics();
      if (!utils::positive(cent_kin_weight))
        return 0.;
      const auto cent_me = computeFactorisedMatrixElement();
      if (!utils::positive(cent_me))
        return 0.;  // avoid computing the fluxes if the matrix element is already null or invalid
      const auto fluxes_weight = part_psgen_->fluxes();
      if (!utils::positive(fluxes_weight))
        return 0.;
      return fluxes_weight * cent_kin_weight * cent_me;
    }

    void FactorisedProcess::fillKinematics() {
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

    //----- utilities

    double FactorisedProcess::that() const {
      //FIXME only works for 2-to-4
      const double that1 = (q1() - pc(0)).mass2();
      const double that2 = (q2() - pc(1)).mass2();
      return 0.5 * (that1 + that2);
    }

    double FactorisedProcess::uhat() const {
      //FIXME only works for 2-to-4
      const double uhat1 = (q1() - pc(1)).mass2();
      const double uhat2 = (q2() - pc(0)).mass2();
      return 0.5 * (uhat1 + uhat2);
    }

    ParametersDescription FactorisedProcess::description() {
      auto desc = Process::description();
      desc.setDescription("Unnamed factorised process");
      desc.add<bool>("storeAlphas", false)
          .setDescription("store the electromagnetic and strong coupling constants to the event content?");
      desc.add<ParametersDescription>("kinematicsGenerator", ParametersDescription().setName<std::string>("2to4"));
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
