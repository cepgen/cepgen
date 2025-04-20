/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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
#include "CepGen/Core/GeneratorWorker.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

using namespace cepgen;

GeneratorWorker::GeneratorWorker(const ParametersList& params) : SteeredObject(params) {}

void GeneratorWorker::setRunParameters(const RunParameters* params) {
  run_params_ = params;
  integrand_ = std::make_unique<ProcessIntegrand>(params);
  CG_DEBUG("GeneratorWorker") << "New generator worker initialised for integration/event generation.\n\t"
                              << "Run parameters at " << dynamic_cast<const void*>(run_params_) << ".";
}

GeneratorWorker::~GeneratorWorker() {
  CG_DEBUG("GeneratorWorker") << "Generator worker destructed. Releasing the parameters at "
                              << dynamic_cast<const void*>(run_params_) << ".";
}

void GeneratorWorker::setIntegrator(const Integrator* integrator) {
  integrator_ = integrator;
  CG_DEBUG("GeneratorWorker:integrator") << "Dim-" << integrand_->size() << " " << integrator_->name()
                                         << " integrator set.";
}

void GeneratorWorker::generate(size_t num_events, const std::function<void(const proc::Process&)>& callback) {
  if (!run_params_)
    throw CG_FATAL("GeneratorWorker:generate") << "No steering parameters specified!";
  callback_proc_ = callback;
  while (run_params_->numGeneratedEvents() < num_events)
    next();
}

bool GeneratorWorker::storeEvent() const {
  CG_TICKER(const_cast<RunParameters*>(run_params_)->timeKeeper());

  if (!integrand_->process().hasEvent())
    return true;

  const auto& event = integrand_->process().event();
  if (const auto num_events_generated = run_params_->numGeneratedEvents();
      (num_events_generated + 1) % run_params_->generation().printEvery() == 0)
    CG_DEBUG("GeneratorWorker:store") << utils::s("event", num_events_generated + 1, true) << " generated.";
  if (callback_proc_)
    callback_proc_(integrand_->process());
  for (const auto& event_exporter : run_params_->eventExportersSequence())
    if (!(*event_exporter << event))
      return false;
  const_cast<RunParameters*>(run_params_)->addGenerationTime(event.metadata("time:total"));
  return true;
}

ParametersDescription GeneratorWorker::description() {
  auto desc = ParametersDescription();
  desc.setDescription("Unnamed generator worker");
  return desc;
}
