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

#include <chrono>

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/GeneratorWorker.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/GeneratorWorkerFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

using namespace cepgen;

Generator::Generator(bool safe_mode) : parameters_(new RunParameters) {
  static bool kInitialised = false;
  if (!kInitialised) {
    cepgen::initialise(safe_mode);
    kInitialised = true;
    CG_DEBUG("Generator:init") << "Generator initialised";
  }
}

Generator::Generator(RunParameters* run_parameters) : parameters_(run_parameters) {}

Generator::~Generator() {
  if (parameters_->timeKeeper())
    CG_INFO("Generator:destructor") << parameters_->timeKeeper()->summary();
}

void Generator::clearRun() {
  CG_DEBUG("Generator:clearRun") << "Run is set to be cleared.";
  worker_ = GeneratorWorkerFactory::get().build(parameters_->generation().parameters().get<ParametersList>("worker"));
  CG_DEBUG("Generator:clearRun") << "Initialised a generator worker with parameters: " << worker_->parameters() << ".";
  // destroy and recreate the integrator instance
  resetIntegrator();

  worker_->setRunParameters(parameters_.get());
  worker_->setIntegrator(integrator_.get());
  cross_section_ = Value{};
  parameters_->prepareRun();
  initialised_ = false;
}

void Generator::parseRunParameters(const std::string& filename) {
  setRunParameters(CardsHandlerFactory::get().buildFromFilename(filename)->parseFile(filename).runParameters());
}

const RunParameters& Generator::runParameters() const {
  if (!parameters_)
    throw CG_FATAL("Generator:runParameters") << "Run parameters object is not yet initialised.";
  return *parameters_;
}

RunParameters& Generator::runParameters() {
  if (!parameters_)
    throw CG_FATAL("Generator:runParameters") << "Run parameters object is not yet initialised.";
  return *parameters_;
}

void Generator::setRunParameters(std::unique_ptr<RunParameters>& run_parameters) {
  parameters_ = std::move(run_parameters);
}

double Generator::computePoint(const std::vector<double>& coordinates) {
  if (!worker_)
    clearRun();
  if (!parameters_->hasProcess())
    throw CG_FATAL("Generator:computePoint") << "Trying to compute a point with no process specified!";
  const size_t ndim = worker_->integrand().process().ndim();
  if (coordinates.size() != ndim)
    throw CG_FATAL("Generator:computePoint")
        << "Invalid phase space dimension (ndim=" << ndim << ", given=" << coordinates.size() << ").";
  const auto res = worker_->integrand().eval(coordinates);
  CG_DEBUG("Generator:computePoint") << "Result for x[" << ndim << "] = " << coordinates << ":\n\t" << res << ".";
  return res;
}

void Generator::computeXsection(double& cross_section, double& err) {
  const auto xsec = computeXsection();
  cross_section = xsec;
  err = xsec.uncertainty();
}

Value Generator::computeXsection() {
  CG_INFO("Generator") << "Starting the computation of the process cross-section.";

  integrate();  // run is cleared here

  if (cross_section_ < 1.e-2)
    CG_INFO("Generator") << "Total cross section: " << cross_section_ * 1.e3 << " fb.";
  else if (cross_section_ < 0.5e3)
    CG_INFO("Generator") << "Total cross section: " << cross_section_ << " pb.";
  else if (cross_section_ < 0.5e6)
    CG_INFO("Generator") << "Total cross section: " << cross_section_ * 1.e-3 << " nb.";
  else if (cross_section_ < 0.5e9)
    CG_INFO("Generator") << "Total cross section: " << cross_section_ * 1.e-6 << " µb.";
  else
    CG_INFO("Generator") << "Total cross section: " << cross_section_ * 1.e-9 << " mb.";

  return cross_section_;
}

void Generator::resetIntegrator() {
  CG_TICKER(parameters_->timeKeeper());
  setIntegrator(IntegratorFactory::get().build(
      parameters_->integrator()));  // create a spec-defined integrator in the current scope
}

void Generator::setIntegrator(std::unique_ptr<Integrator> integrator) {
  CG_TICKER(parameters_->timeKeeper());
  if (integrator)  // copy the integrator instance (or create it if unspecified) in the current scope
    integrator_ = std::move(integrator);
  else
    resetIntegrator();
  CG_INFO("Generator:integrator") << "Generator will use a " << integrator_->name() << "-type integrator.";
}

Integrator& Generator::integrator() const {
  if (!integrator_)
    throw CG_FATAL("Generator:integrator") << "Uninitialised integrator object.";
  return *integrator_;
}

void Generator::integrate() {
  CG_TICKER(parameters_->timeKeeper());

  clearRun();
  if (!parameters_->hasProcess())
    throw CG_FATAL("Generator:integrate") << "Trying to integrate while no process is specified!";
  const size_t num_dimensions = worker_->integrand().process().ndim();
  if (num_dimensions == 0)
    throw CG_FATAL("Generator:integrate") << "Invalid phase space dimension. "
                                          << "At least one integration variable is required!";

  CG_DEBUG("Generator:integrate") << "New integrator instance created for " << num_dimensions
                                  << "-dimensional integration.";

  if (!integrator_)
    throw CG_FATAL("Generator:integrate") << "No integrator object was declared for the generator!";

  cross_section_ = integrator_->integrate(worker_->integrand());

  CG_DEBUG("Generator:integrate") << "Computed cross section: (" << cross_section_ << ") pb.";

  // now that the cross-section has been computed, feed it to the event modification algorithms...
  for (const auto& event_modifier : parameters_->eventModifiersSequence())
    event_modifier->setCrossSection(cross_section_);
  // ...and to the event storage algorithms
  for (const auto& event_exporter : parameters_->eventExportersSequence())
    event_exporter->setCrossSection(cross_section_);
}

void Generator::initialise() {
  if (!parameters_)
    throw CG_FATAL("Generator:generate") << "No steering parameters specified!";

  CG_TICKER(parameters_->timeKeeper());

  if (!worker_)  // if no worker is found, launch a new integration/run preparation
    integrate();

  // prepare the run parameters for event generation
  parameters_->initialiseModules();
  worker_->initialise();
  initialised_ = true;
}

const Event& Generator::next() {
  if (!worker_ || !initialised_)
    initialise();
  size_t num_try = 0;
  while (!worker_->next()) {
    if (num_try++ > 5)
      throw CG_FATAL("Generator:next") << "Failed to generate the next event!";
  }
  return worker_->integrand().process().event();
}

void Generator::generate(size_t num_events, const std::function<void(const proc::Process&)>& callback) {
  CG_TICKER(parameters_->timeKeeper());

  if (!worker_ || !initialised_)
    initialise();

  // if invalid argument, retrieve from runtime parameters
  if (num_events < 1) {
    if (parameters_->generation().targetLuminosity() > 0.) {
      num_events = std::ceil(parameters_->generation().targetLuminosity() * cross_section_);
      CG_INFO("Generator") << "Target luminosity: " << parameters_->generation().targetLuminosity() << " pb-1.";
    } else
      num_events = parameters_->generation().maxGen();
  }

  CG_INFO("Generator") << utils::s("event", num_events, true) << " will be generated.";

  const utils::Timer tmr;

  worker_->generate(num_events, callback);  // launch the event generation

  const double generation_time = tmr.elapsed();
  const double rate_ms = (parameters_->numGeneratedEvents() > 0)
                             ? generation_time * 1.e3 /* s -> ms */ / parameters_->numGeneratedEvents()
                             : 0.;
  const double equivalent_luminosity = parameters_->numGeneratedEvents() / crossSection();
  CG_INFO("Generator") << utils::s("event", parameters_->numGeneratedEvents()) << " generated in " << generation_time
                       << " s "
                       << "(" << rate_ms << " ms/event).\n\t"
                       << "Equivalent luminosity: " << utils::format("%g", equivalent_luminosity) << " pb^-1.";
}

void Generator::generate(size_t num_events, const std::function<void(const Event&, size_t)>& callback) {
  generate(num_events, [this, &callback](const proc::Process& process) {
    if (callback)
      callback(process.event(), parameters_->numGeneratedEvents());
  });
}
