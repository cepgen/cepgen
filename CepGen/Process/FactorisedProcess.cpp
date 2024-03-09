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
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PhaseSpaceGenerator.h"
#include "CepGen/Utils/Math.h"

namespace cepgen {
  namespace proc {
    FactorisedProcess::FactorisedProcess(const ParametersList& params, const spdgids_t& central)
        : Process(params),
          psgen_(PhaseSpaceGeneratorFactory::get().build(
              steer<ParametersList>("kinematicsGenerator")
                  .set("ids", std::vector<int>(central.begin(), central.end())))),
          store_alphas_(steer<bool>("storeAlphas")) {
      event().map()[Particle::CentralSystem].resize(central.size());
    }

    FactorisedProcess::FactorisedProcess(const FactorisedProcess& proc)
        : Process(proc),
          psgen_(PhaseSpaceGeneratorFactory::get().build(proc.psgen_->parameters())),
          store_alphas_(proc.store_alphas_) {}

    void FactorisedProcess::addEventContent() {
      CG_ASSERT(psgen_);
      const auto cent_pdgids = psgen_->central();
      Process::setEventContent({{Particle::IncomingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
                                {Particle::IncomingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
                                {Particle::OutgoingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
                                {Particle::OutgoingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
                                {Particle::CentralSystem, spdgids_t(cent_pdgids.begin(), cent_pdgids.end())}});
    }

    void FactorisedProcess::prepareKinematics() {
      if (!psgen_)
        throw CG_FATAL("FactorisedProcess:prepareKinematics")
            << "Phase space generator not set. Please check your process initialisation procedure, as you might "
               "be doing something irregular.";
      psgen_->initialise(this);

      event().oneWithRole(Particle::Parton1).setPdgId(psgen_->partons().at(0));
      event().oneWithRole(Particle::Parton2).setPdgId(psgen_->partons().at(1));

      CG_DEBUG("FactorisedProcess:prepareKinematics") << "Partons: " << psgen_->partons() << ", "
                                                      << "central system: " << psgen_->central() << ". " << event();

      // register all process-dependent variables
      prepareFactorisedPhaseSpace();

      // register the outgoing remnants' variables
      if (!kinematics().incomingBeams().positive().elastic())
        defineVariable(mX2(), Mapping::square, kinematics().cuts().remnants.mx, "Positive-z beam remnant squared mass");
      if (!kinematics().incomingBeams().negative().elastic())
        defineVariable(mY2(), Mapping::square, kinematics().cuts().remnants.mx, "Negative-z beam remnant squared mass");
    }

    double FactorisedProcess::computeWeight() {
      if (!psgen_->generate())
        return 0.;
      if (const auto cent_weight = computeFactorisedMatrixElement(); utils::positive(cent_weight))
        return cent_weight * psgen_->weight();
      return 0.;
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

    double FactorisedProcess::that() const { return psgen_->that(); }

    double FactorisedProcess::uhat() const { return psgen_->uhat(); }

    ParametersDescription FactorisedProcess::description() {
      auto desc = Process::description();
      desc.setDescription("Unnamed factorised process");
      desc.add<bool>("storeAlphas", false)
          .setDescription("store the electromagnetic and strong coupling constants to the event content?");
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
