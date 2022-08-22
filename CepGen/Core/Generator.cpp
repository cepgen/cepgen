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

#include <chrono>

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/GeneratorWorker.h"
#include "CepGen/Generator.h"
#include "CepGen/Integration/GridParameters.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

namespace cepgen {
  Generator::Generator(bool safe_mode) : parameters_(new Parameters) {
    static bool init = false;
    if (!init) {
      initialise(safe_mode);
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
    worker_.reset(new GeneratorWorker(const_cast<const Parameters*>(parameters_.get())));
    // destroy and recreate the integrator instance
    setIntegrator(nullptr);
    worker_->setIntegrator(integrator_.get());
    result_ = result_error_ = -1.;
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
      throw CG_FATAL("Generator:computePoint") << "Invalid phase space dimension (ndim=" << ndim << ")!";
    double res = worker_->integrand().eval(coord);
    CG_DEBUG("Generator:computePoint") << "Result for x[" << ndim << "] = " << coord << ":\n\t" << res << ".";
    return res;
  }

  void Generator::computeXsection(double& cross_section, double& err) {
    CG_INFO("Generator") << "Starting the computation of the process cross-section.";

    integrate();  // run is cleared here

    cross_section = result_;
    err = result_error_;

    if (cross_section < 1.e-2)
      CG_INFO("Generator") << "Total cross section: " << cross_section * 1.e3 << " +/- " << err * 1.e3 << " fb.";
    else if (cross_section < 0.5e3)
      CG_INFO("Generator") << "Total cross section: " << cross_section << " +/- " << err << " pb.";
    else if (cross_section < 0.5e6)
      CG_INFO("Generator") << "Total cross section: " << cross_section * 1.e-3 << " +/- " << err * 1.e-3 << " nb.";
    else if (cross_section < 0.5e9)
      CG_INFO("Generator") << "Total cross section: " << cross_section * 1.e-6 << " +/- " << err * 1.e-6 << " µb.";
    else
      CG_INFO("Generator") << "Total cross section: " << cross_section * 1.e-9 << " +/- " << err * 1.e-9 << " mb.";
  }

  void Generator::setIntegrator(std::unique_ptr<Integrator> integ) {
    CG_TICKER(parameters_->timeKeeper());

    if (!integ) {
      if (parameters_->par_integrator.name<std::string>().empty())
        parameters_->par_integrator.setName<std::string>("Vegas");
      integ = IntegratorFactory::get().build(parameters_->par_integrator);
    }
    integrator_ = std::move(integ);
    integrator_->setIntegrand(worker_->integrand());
    if (worker_)
      worker_->setIntegrator(integrator_.get());
    CG_INFO("Generator:integrator") << "Generator will use a " << integrator_->name() << "-type integrator.";
  }

  void Generator::integrate() {
    CG_TICKER(parameters_->timeKeeper());

    clearRun();
    if (!parameters_->hasProcess())
      throw CG_FATAL("Generator:integrate") << "Trying to integrate while no process is specified!";
    const size_t ndim = worker_->integrand().process().ndim();
    if (ndim == 0)
      throw CG_FATAL("Generator:computePoint") << "Invalid phase space dimension. "
                                               << "At least one integration variable is required!";

    CG_DEBUG("Generator:integrate") << "New integrator instance created for " << ndim << "-dimensional integration.";

    integrator_->integrate(result_, result_error_);

    CG_DEBUG("Generator:integrate") << "Computed cross section: (" << result_ << " +- " << result_error_ << ") pb.";

    // now that the cross section has been computed, feed it to the event modification algorithms...
    for (auto& mod : parameters_->eventModifiersSequence())
      mod->setCrossSection(result_, result_error_);
    // ...and to the event storage algorithms
    for (auto& mod : parameters_->outputModulesSequence())
      mod->setCrossSection(result_, result_error_);
  }

  const Event& Generator::generateOneEvent(Event::callback callback) { return next(callback); }

  void Generator::initialise() {
    CG_TICKER(parameters_->timeKeeper());

    if (!parameters_)
      throw CG_FATAL("Generator:generate") << "No steering parameters specified!";
    if (!worker_)
      integrate();

    for (auto& mod : parameters_->outputModulesSequence())
      mod->initialise(*parameters_);

    initialised_ = true;
  }

  const Event& Generator::next(Event::callback callback) {
    if (!worker_ || !initialised_)
      initialise();
    return worker_->next(callback);
  }

  void Generator::generate(size_t num_events, Event::callback callback) {
    CG_TICKER(parameters_->timeKeeper());

    if (!initialised_)
      initialise();

    //--- if invalid argument, retrieve from runtime parameters
    if (num_events < 1) {
      if (parameters_->generation().targetLuminosity() > 0.) {
        num_events = std::ceil(parameters_->generation().targetLuminosity() * result_);
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
}  // namespace cepgen
