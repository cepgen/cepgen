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
    FactorisedProcess::FactorisedProcess(const ParametersList& params, const pdgids_t& central)
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
      Process::setEventContent({{Particle::IncomingBeam1, {kinematics().incomingBeams().positive().pdgId()}},
                                {Particle::IncomingBeam2, {kinematics().incomingBeams().negative().pdgId()}},
                                {Particle::OutgoingBeam1, {kinematics().incomingBeams().positive().pdgId()}},
                                {Particle::OutgoingBeam2, {kinematics().incomingBeams().negative().pdgId()}},
                                {Particle::CentralSystem, psgen_->central()}});
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
      if (const auto ps_weight = psgen_->generate(); utils::positive(ps_weight))
        return ps_weight * computeFactorisedMatrixElement();
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
      return desc;
    }
  }  // namespace proc
}  // namespace cepgen
