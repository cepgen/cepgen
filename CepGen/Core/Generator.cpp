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

#include <chrono>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/GeneratorWorker.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Modules/GeneratorWorkerFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

namespace cepgen {
  Generator::Generator(bool safe_mode) : parameters_(new Parameters) {
    static bool init = false;
    if (!init) {
      cepgen::initialise(safe_mode);
      init = true;
      CG_DEBUG("Generator:init") << "Generator initialised";
    }
    //--- random number initialization
    std::chrono::system_clock::time_point time = std::chrono::system_clock::now();
    srandom(time.time_since_epoch().count());
  }

  Generator::Generator(Parameters* ip) : parameters_(ip) {}

  Generator::~Generator() {
    if (parameters_->timeKeeper())
      CG_INFO("Generator:destructor") << parameters_->timeKeeper()->summary();
  }

  void Generator::clearRun() {
    CG_DEBUG("Generator:clearRun") << "Run is set to be cleared.";
    worker_ = GeneratorWorkerFactory::get().build(parameters_->generation().parameters().get<ParametersList>("worker"));
    CG_DEBUG("Generator:clearRun") << "Initialised a generator worker with parameters: " << worker_->parameters()
                                   << ".";
    // destroy and recreate the integrator instance
    if (!integrator_)
      resetIntegrator();

    worker_->setRuntimeParameters(const_cast<const Parameters*>(parameters_.get()));
    worker_->setIntegrator(integrator_.get());
    xsect_ = Value{-1., -1.};
    parameters_->prepareRun();
  }

  Parameters& Generator::parametersRef() { return *parameters_; }

  void Generator::setParameters(Parameters* ip) { parameters_.reset(ip); }

  double Generator::computePoint(const std::vector<double>& coord) {
    if (!worker_)
      clearRun();
    if (!parameters_->hasProcess())
      throw CG_FATAL("Generator:computePoint") << "Trying to compute a point with no process specified!";
    const size_t ndim = worker_->integrand().process().ndim();
    if (coord.size() != ndim)
      throw CG_FATAL("Generator:computePoint")
          << "Invalid phase space dimension (ndim=" << ndim << ", given=" << coord.size() << ").";
    double res = worker_->integrand().eval(coord);
    CG_DEBUG("Generator:computePoint") << "Result for x[" << ndim << "] = " << coord << ":\n\t" << res << ".";
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

    if (xsect_ < 1.e-2)
      CG_INFO("Generator") << "Total cross section: " << xsect_ * 1.e3 << " fb.";
    else if (xsect_ < 0.5e3)
      CG_INFO("Generator") << "Total cross section: " << xsect_ << " pb.";
    else if (xsect_ < 0.5e6)
      CG_INFO("Generator") << "Total cross section: " << xsect_ * 1.e-3 << " nb.";
    else if (xsect_ < 0.5e9)
      CG_INFO("Generator") << "Total cross section: " << xsect_ * 1.e-6 << " µb.";
    else
      CG_INFO("Generator") << "Total cross section: " << xsect_ * 1.e-9 << " mb.";

    return xsect_;
  }

  void Generator::resetIntegrator() {
    CG_TICKER(parameters_->timeKeeper());
    // create a spec-defined integrator in the current scope
    setIntegrator(IntegratorFactory::get().build(parameters_->par_integrator));
  }

  void Generator::setIntegrator(std::unique_ptr<Integrator> integ) {
    CG_TICKER(parameters_->timeKeeper());
    // copy the integrator instance (or create it if unspecified) in the current scope
    if (!integ)
      resetIntegrator();
    else
      integrator_ = std::move(integ);
    CG_INFO("Generator:integrator") << "Generator will use a " << integrator_->name() << "-type integrator.";
  }

  void Generator::integrate() {
    CG_TICKER(parameters_->timeKeeper());

    clearRun();

    if (!parameters_->hasProcess())
      throw CG_FATAL("Generator:integrate") << "Trying to integrate while no process is specified!";
    const size_t ndim = worker_->integrand().process().ndim();
    if (ndim == 0)
      throw CG_FATAL("Generator:integrate") << "Invalid phase space dimension. "
                                            << "At least one integration variable is required!";

    CG_DEBUG("Generator:integrate") << "New integrator instance created for " << ndim << "-dimensional integration.";

    if (!integrator_)
      throw CG_FATAL("Generator:integrate") << "No integrator object was declared for the generator!";

    xsect_ = integrator_->integrate(worker_->integrand());

    CG_DEBUG("Generator:integrate") << "Computed cross section: (" << xsect_ << ") pb.";

    // now that the cross section has been computed, feed it to the event modification algorithms...
    for (auto& mod : parameters_->eventModifiersSequence())
      mod->setCrossSection(xsect_);
    // ...and to the event storage algorithms
    for (auto& mod : parameters_->eventExportersSequence())
      mod->setCrossSection(xsect_);
  }

  void Generator::initialise() {
    if (initialised_)
      return;
    if (!parameters_)
      throw CG_FATAL("Generator:generate") << "No steering parameters specified!";

    CG_TICKER(parameters_->timeKeeper());

    // if no worker is found, launch a new integration/run preparation
    if (!worker_)
      integrate();

    // prepare the run parameters for event generation
    parameters_->initialise();
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

    //--- if invalid argument, retrieve from runtime parameters
    if (num_events < 1) {
      if (parameters_->generation().targetLuminosity() > 0.) {
        num_events = std::ceil(parameters_->generation().targetLuminosity() * xsect_);
        CG_INFO("Generator") << "Target luminosity: " << parameters_->generation().targetLuminosity() << " pb-1.";
      } else
        num_events = parameters_->generation().maxGen();
    }

    CG_INFO("Generator") << utils::s("event", num_events, true) << " will be generated.";

    const utils::Timer tmr;

    //--- launch the event generation

    worker_->generate(num_events, callback);

    const double gen_time_s = tmr.elapsed();
    const double rate_ms =
        (parameters_->numGeneratedEvents() > 0) ? gen_time_s / parameters_->numGeneratedEvents() * 1.e3 : 0.;
    const double equiv_lumi = parameters_->numGeneratedEvents() / crossSection();
    CG_INFO("Generator") << utils::s("event", parameters_->numGeneratedEvents()) << " generated in " << gen_time_s
                         << " s "
                         << "(" << rate_ms << " ms/event).\n\t"
                         << "Equivalent luminosity: " << utils::format("%g", equiv_lumi) << " pb^-1.";
  }

  void Generator::generate(size_t num_events, const std::function<void(const Event&, size_t)>& callback) {
    generate(num_events, [&](const proc::Process& proc) {
      if (callback)
        callback(proc.event(), parameters_->numGeneratedEvents());
    });
  }
}  // namespace cepgen
