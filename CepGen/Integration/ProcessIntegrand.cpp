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

#include <numeric>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
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
    if (!process_->hasEvent())
      return weight;

    //--- prepare the event content
    auto* event = process_->eventPtr();

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
                               [&bws, &event](double init, const auto& tam) {
                                 return init * (*tam)(bws.get(*event, tam->variables().at(0)));
                               });
    } catch (const Exception&) {
      throw CG_FATAL("ProcessIntegrand") << "Failed to apply taming function(s) taming!";
    }

    if (weight <= 0.)
      return 0.;

    //--- set the CepGen part of the event generation
    if (storage_)
      event->time_generation = (float)tmr_->elapsed();

    //--- trigger all event modification algorithms
    if (!params_->eventModifiersSequence().empty()) {
      double br = -1.;
      for (auto& mod : params_->eventModifiersSequence()) {
        if (!mod->run(*event, br, storage_) || br == 0.)
          return 0.;
        weight *= br;  // branching fraction for all decays
      }
    }

    //--- apply cuts on final state system (after hadronisation!)
    //    (polish your cuts, as this might be very time-consuming...)

    if (!params_->kinematics().cuts().central.contain((*event)(Particle::CentralSystem)))
      return 0.;
    if (!params_->kinematics().cuts().central_particles.empty())
      for (const auto& part : (*event)(Particle::CentralSystem)) {
        // retrieve all cuts associated to this final state particle in the central system
        if (params_->kinematics().cuts().central_particles.count(part.pdgId()) > 0 &&
            !params_->kinematics().cuts().central_particles.at(part.pdgId()).contain({part}))
          return 0.;
      }
    if (!params_->kinematics().cuts().remnants.contain((*event)(Particle::OutgoingBeam1), event))
      return 0.;
    if (!params_->kinematics().cuts().remnants.contain((*event)(Particle::OutgoingBeam2), event))
      return 0.;

    //--- store the last event in parameters block for a later usage
    if (storage_) {
      event->weight = (float)weight;
      event->time_total = (float)tmr_->elapsed();

      CG_DEBUG_LOOP("ProcessIntegrand") << "[process " << std::hex << (void*)process_.get() << std::dec << "]\n\t"
                                        << "Generation time: " << event->time_generation * 1.e3 << " ms\n\t"
                                        << "Total time (gen+hadr+cuts): " << event->time_total * 1.e3 << " ms";
    }

    //--- a bit of useful debugging
    CG_DEBUG_LOOP("ProcessIntegrand") << "f value for dim-" << x.size() << " point " << x << ": " << weight << ".";

    return weight;
  }
}  // namespace cepgen
