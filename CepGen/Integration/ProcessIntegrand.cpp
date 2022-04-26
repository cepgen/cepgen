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

#include <numeric>

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/EventBrowser.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/TimeKeeper.h"

namespace cepgen {
  ProcessIntegrand::ProcessIntegrand(const Parameters* params) : params_(params), tmr_(new utils::Timer) {
    if (!params_) {
      CG_WARNING("ProcessIntegrand") << "Invalid runtime parameters specified.";
      return;
    }
    if (!params_->hasProcess()) {
      CG_WARNING("ProcessIntegrand") << "No process defined in runtime parameters.";
      return;
    }
    //--- each integrand object has its own clone of the process
    process_ = params_->process().clone();
    CG_DEBUG("ProcessIntegrand") << "Process " << process_->name() << " successfully cloned from base process "
                                 << params_->process().name() << ".";

    //--- process-specific phase space definition
    process_->prepareKinematics();
    CG_DEBUG("ProcessIntegrand").log([this](auto& log) { process_->dumpVariables(&log.stream()); });

    //--- prepare the event content
    if (process_->hasEvent())
      event_ = &process_->event();

    //--- first-run preparation
    CG_DEBUG("ProcessIntegrand") << "Preparing all variables for a new run.";
    const auto& kin = process_->kinematics();
    CG_DEBUG("ProcessIntegrand").log([&](auto& dbg) {
      dbg << "Run started for " << process_->name() << " process " << std::hex << (void*)process_.get() << std::dec
          << ".\n\t";
      const auto& beams = kin.incomingBeams();
      dbg << "Process mode considered: " << beams.mode() << "\n\t"
          << "  positive-z beam: " << beams.positive() << "\n\t"
          << "  negative-z beam: " << beams.negative();
      if (beams.structureFunctions())
        dbg << "\n\t  structure functions: " << beams.structureFunctions();
    });
    if (process_->hasEvent())
      process_->clearEvent();

    CG_DEBUG("ProcessIntegrand") << "New integrand object defined for process \"" << process_->name() << "\".";
  }

  ProcessIntegrand::~ProcessIntegrand() { CG_DEBUG("ProcessIntegrand") << "Destructor called"; }

  size_t ProcessIntegrand::size() const {
    if (!process_)
      throw CG_FATAL("ProcessIntegrand:size") << "Process was not properly cloned!";
    return process_->ndim();
  }

  double ProcessIntegrand::eval(const std::vector<double>& x) {
    CG_TICKER(const_cast<Parameters*>(params_)->timeKeeper());

    //--- start the timer
    tmr_->reset();

    //--- specify the phase space point to probe and calculate weight
    double weight = process_->weight(x);

    //--- invalidate any unphysical behaviour
    if (weight <= 0.)
      return 0.;

    //--- speed up the integration process if no event is to be generated
    if (!event_)
      return weight;
    if (!storage_ && !params_->eventModifiersSequence().empty() && !params_->tamingFunctions().empty() &&
        params_->kinematics().cuts().central_particles.empty())
      return weight;

    //--- fill in the process' Event object
    process_->fillKinematics();

    //--- once the kinematics variables have been populated, can apply the
    //    collection of taming functions
    try {
      utils::EventBrowser bws;
      weight = std::accumulate(params_->tamingFunctions().begin(),
                               params_->tamingFunctions().end(),
                               weight,
                               [&bws, this](double init, const auto& tam) {
                                 return init * (*tam)(bws.get(*event_, tam->variables().at(0)));
                               });
    } catch (const Exception&) {
      throw CG_FATAL("ProcessIntegrand") << "Failed to apply taming function(s) taming!";
    }

    if (weight <= 0.)
      return 0.;

    //--- set the CepGen part of the event generation
    if (storage_)
      event_->time_generation = (float)tmr_->elapsed();

    //--- trigger all event modification algorithms
    if (!params_->eventModifiersSequence().empty()) {
      double br = -1.;
      for (auto& mod : params_->eventModifiersSequence()) {
        if (!mod->run(*event_, br, storage_) || br == 0.)
          return 0.;
        weight *= br;  // branching fraction for all decays
      }
    }

    //--- apply cuts on final state system (after hadronisation!)
    //    (polish your cuts, as this might be very time-consuming...)

    if (!params_->kinematics().cuts().central_particles.empty())
      for (const auto& part : (*event_)(Particle::CentralSystem)) {
        // retrieve all cuts associated to this final state particle in the
        // central system
        if (params_->kinematics().cuts().central_particles.count(part.pdgId()) == 0)
          continue;
        const auto& cuts_pdgid = params_->kinematics().cuts().central_particles.at(part.pdgId());
        // apply these cuts on the given particle
        if (!cuts_pdgid.pt_single().contains(part.momentum().pt()))
          return 0.;
        if (!cuts_pdgid.energy_single().contains(part.momentum().energy()))
          return 0.;
        if (!cuts_pdgid.eta_single().contains(part.momentum().eta()))
          return 0.;
        if (!cuts_pdgid.rapidity_single().contains(part.momentum().rapidity()))
          return 0.;
      }
    const auto& remn_cut = params_->kinematics().cuts().remnants;
    for (const auto& system : {Particle::OutgoingBeam1, Particle::OutgoingBeam2})
      for (const auto& part : (*event_)(system)) {
        if (part.status() != Particle::Status::FinalState)
          continue;
        if (!remn_cut.xi().contains(1. - part.momentum().pz() / (*event_)[*part.mothers().begin()].momentum().pz()))
          return 0.;
        if (!remn_cut.yj().contains(fabs(part.momentum().rapidity())))
          return 0.;
      }

    //--- store the last event in parameters block for a later usage
    if (storage_) {
      event_->weight = (float)weight;
      event_->time_total = (float)tmr_->elapsed();

      CG_DEBUG_LOOP("ProcessIntegrand") << "[process " << std::hex << (void*)process_.get() << std::dec << "]\n\t"
                                        << "Generation time: " << event_->time_generation * 1.e3 << " ms\n\t"
                                        << "Total time (gen+hadr+cuts): " << event_->time_total * 1.e3 << " ms";
    }

    //--- a bit of useful debugging
    CG_DEBUG_LOOP("ProcessIntegrand") << "f value for dim-" << x.size() << " point " << x << ": " << weight << ".";

    return weight;
  }
}  // namespace cepgen
