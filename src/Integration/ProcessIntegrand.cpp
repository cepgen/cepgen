/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/TimeKeeper.h"

using namespace cepgen;

ProcessIntegrand::ProcessIntegrand(const proc::Process& process)
    : run_parameters_(new RunParameters), timer_(new utils::Timer) {
  setProcess(process);
}

ProcessIntegrand::ProcessIntegrand(const RunParameters* run_parameters)
    : run_parameters_(run_parameters), timer_(new utils::Timer) {
  if (!run_parameters_)
    throw CG_FATAL("ProcessIntegrand") << "Invalid runtime parameters specified.";
  if (!run_parameters_->hasProcess())
    throw CG_FATAL("ProcessIntegrand") << "No process defined in runtime parameters.";
  setProcess(run_parameters_->process());
}

size_t ProcessIntegrand::size() const { return process().ndim(); }

void ProcessIntegrand::setProcess(const proc::Process& original_process) {
  process_ = original_process.clone();  // each integrand object has its own clone of the process
  // NOTE: kinematics is already set by the process copy constructor
  process().kinematics().setParameters(
      original_process.kinematics().parameters());  // override default kinematics with the one defined in mother process
  CG_DEBUG("ProcessIntegrand:setProcess")
      << "New '" << process().name() << "' process cloned from '" << original_process.name()
      << "' process. New kinematics: " << process().kinematics().parameters() << ".";

  // first-run preparation
  CG_DEBUG("ProcessIntegrand:setProcess").log([this](auto& log) {
    log << "Run started for " << process().name() << " process " << std::hex << dynamic_cast<void*>(process_.get())
        << std::dec << ".\n\t";
    const auto& beams = process().kinematics().incomingBeams();
    log << "Process mode considered: " << beams.mode() << "\n\t"
        << "  positive-z beam: " << beams.positive() << "\n\t"
        << "  negative-z beam: " << beams.negative();
    if (!beams.structureFunctions().empty())
      log << "\n\t  structure functions: " << beams.structureFunctions();
    process().dumpVariables(&log.stream());
  });
  process().initialise();

  CG_DEBUG("ProcessIntegrand:setProcess")
      << "Process integrand defined for dimension-" << size() << " process '" << process().name() << "'.";
}

proc::Process& ProcessIntegrand::process() {
  if (!process_)
    throw CG_FATAL("ProcessIntegrand:process") << "Process was not properly cloned!";
  return *process_;
}

const proc::Process& ProcessIntegrand::process() const {
  if (!process_)
    throw CG_FATAL("ProcessIntegrand:process") << "Process was not properly cloned!";
  return *process_;
}

double ProcessIntegrand::eval(const std::vector<double>& x) {
  CG_TICKER(const_cast<RunParameters*>(run_parameters_)->timeKeeper());
  timer_->reset();  // start the timer

  process().clearEvent();
  auto weight = process().weight(x);  // specify the phase space point to probe and calculate weight
  if (!utils::positive(weight))       // invalidate any unphysical behaviour
    return 0.;

  if (!process_->hasEvent())  // speed up the integration process if no event is to be generated
    return weight;
  process_->setKinematics();           // fill in the process' Event object
  auto* event = process_->eventPtr();  // prepare the event content

  // once kinematics variables computed, can apply taming functions
  for (const auto& taming_function : run_parameters_->tamingFunctions())
    if (const auto val = (*taming_function)(bws_.get(*event, taming_function->variables().at(0))) != 0.)
      weight *= val;
    else
      return 0.;

  if (storage_)
    event->metadata["time:generation"] = timer_->elapsed();  // pure CepGen part of the event generation

  {  // run all event modification algorithms
    double branching_ratio = -1.;
    for (auto& event_modifier : run_parameters_->eventModifiersSequence()) {
      if (!event_modifier->run(*event, branching_ratio, !storage_) || branching_ratio == 0.)
        return 0.;
      weight *= branching_ratio;  // branching fraction for all decays
    }
  }
  const auto& kinematics = process_->kinematics();
  if (!kinematics.cuts().central.contain((*event)(Particle::Role::CentralSystem)))
    // apply cuts on final state system (after event modification algorithms)
    // (polish your cuts, as this might be very time-consuming...)
    return 0.;
  for (const auto& part : (*event)(Particle::Role::CentralSystem))
    // retrieve all cuts associated to this final state particle in the central system
    if (kinematics.cuts().central_particles.count(part.pdgId()) > 0 &&
        !kinematics.cuts().central_particles.at(part.pdgId()).contain({part}))
      return 0.;
  if (!kinematics.incomingBeams().positive().elastic() &&
      !kinematics.cuts().remnants.contain((*event)(Particle::Role::OutgoingBeam1), event))
    return 0.;
  if (!kinematics.incomingBeams().negative().elastic() &&
      !kinematics.cuts().remnants.contain((*event)(Particle::Role::OutgoingBeam2), event))
    return 0.;

  if (storage_) {  // add generation metadata to the event
    event->metadata["weight"] = weight;
    event->metadata["time:total"] = timer_->elapsed();
  }

  CG_DEBUG_LOOP("ProcessIntegrand") << "[process " << std::hex << dynamic_cast<void*>(process_.get()) << std::dec
                                    << "]\n\t"
                                    << "functional value for dimension-" << x.size() << " point " << x << ": " << weight
                                    << ".\n\t"
                                    << "Generation time: " << event->metadata("time:generation") * 1.e3 << " ms\n\t"
                                    << "Total time (gen+mod+cuts): " << event->metadata("time:total") * 1.e3 << " ms";
  return weight;
}
