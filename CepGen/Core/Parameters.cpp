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

#include <iomanip>

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

namespace cepgen {
  Parameters::Parameters() {}

  Parameters::Parameters(Parameters& param)
      : par_integrator(param.par_integrator),
        process_(std::move(param.process_)),
        evt_modifiers_(std::move(param.evt_modifiers_)),
        out_modules_(std::move(param.out_modules_)),
        taming_functions_(std::move(param.taming_functions_)),
        total_gen_time_(param.total_gen_time_),
        num_gen_events_(param.num_gen_events_),
        generation_(param.generation_),
        tmr_(std::move(param.tmr_)) {}

  Parameters::Parameters(const Parameters& param)
      : par_integrator(param.par_integrator),
        total_gen_time_(param.total_gen_time_),
        num_gen_events_(param.num_gen_events_),
        generation_(param.generation_) {}

  Parameters::~Parameters()  // required for unique_ptr initialisation!
  {}

  Parameters& Parameters::operator=(Parameters param) {
    par_integrator = param.par_integrator;
    process_ = std::move(param.process_);
    evt_modifiers_ = std::move(param.evt_modifiers_);
    out_modules_ = std::move(param.out_modules_);
    taming_functions_ = std::move(param.taming_functions_);
    total_gen_time_ = param.total_gen_time_;
    num_gen_events_ = param.num_gen_events_;
    generation_ = param.generation_;
    tmr_ = std::move(param.tmr_);
    return *this;
  }

  void Parameters::prepareRun() {
    if (tmr_)
      tmr_->clear();
    CG_TICKER(tmr_.get());

    //--- clear the run statistics
    total_gen_time_ = 0.;
    num_gen_events_ = 0ul;
  }

  void Parameters::setTimeKeeper(utils::TimeKeeper* kpr) { tmr_.reset(kpr); }

  void Parameters::addGenerationTime(double gen_time) {
    total_gen_time_ += gen_time;
    num_gen_events_++;
  }

  proc::Process& Parameters::process() { return *process_.get(); }

  const proc::Process& Parameters::process() const { return *process_.get(); }

  std::string Parameters::processName() const {
    if (!process_)
      return "no process";
    return process_->name();
  }

  void Parameters::clearProcess() { process_.release(); }

  void Parameters::setProcess(std::unique_ptr<proc::Process> proc) { process_ = std::move(proc); }

  void Parameters::setProcess(proc::Process* proc) {
    if (!proc)
      throw CG_FATAL("Parameters") << "Trying to clone an invalid process!";
    process_.reset(proc);
  }

  const Kinematics& Parameters::kinematics() const {
    if (!process_)
      throw CG_FATAL("Parameters") << "Process must be defined before its kinematics is retrieved!";
    return process_->kinematics();
  }

  EventModifier& Parameters::eventModifier(size_t i) { return *evt_modifiers_.at(i); }

  void Parameters::clearEventModifiersSequence() { evt_modifiers_.clear(); }

  void Parameters::addModifier(std::unique_ptr<EventModifier> mod) {
    evt_modifiers_.emplace_back(std::move(mod));
    (*evt_modifiers_.rbegin())->setRuntimeParameters(*this);
  }

  void Parameters::addModifier(EventModifier* mod) {
    std::unique_ptr<EventModifier> modifier(mod);
    modifier->setRuntimeParameters(*this);
    evt_modifiers_.emplace_back(std::move(modifier));
  }

  io::ExportModule& Parameters::outputModule(size_t i) { return *out_modules_.at(i); }

  void Parameters::clearOutputModulesSequence() { out_modules_.clear(); }

  void Parameters::addOutputModule(std::unique_ptr<io::ExportModule> mod) { out_modules_.emplace_back(std::move(mod)); }

  void Parameters::addOutputModule(io::ExportModule* mod) {
    out_modules_.emplace_back(std::unique_ptr<io::ExportModule>(mod));
  }

  void Parameters::addTamingFunction(std::unique_ptr<utils::Functional> fct) {
    taming_functions_.emplace_back(std::move(fct));
  }

  std::ostream& operator<<(std::ostream& os, const Parameters* param) {
    const int wb = 90, wt = 33;

    os << std::left << "\n"
       << std::setfill('_') << std::setw(wb + 3) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill(' ') << "\n\n";
    if (param->process_)
      os << std::setw(wt) << "Process to generate:\n"
         << proc::ProcessFactory::get()
                .describeParameters(param->process().name(), param->process().parameters())
                .describe(1)
         << "\n\n";
    os << std::setw(wt) << "Events generation? " << utils::yesno(param->generation_.enabled()) << "\n"
       << std::setw(wt) << "Number of events to generate" << utils::boldify(param->generation_.maxGen()) << "\n";
    if (param->generation_.numThreads() > 1)
      os << std::setw(wt) << "Number of threads" << param->generation_.numThreads() << "\n";
    os << std::setw(wt) << "Number of points to try per bin" << param->generation_.numPoints() << "\n"
       << std::setw(wt) << "Verbosity level " << utils::Logger::get().level << "\n";
    if (!param->evt_modifiers_.empty() || param->out_modules_.empty() || !param->taming_functions_.empty())
      os << "\n"
         << std::setfill('-') << std::setw(wb + 6) << utils::boldify(" Event treatment ") << std::setfill(' ')
         << "\n\n";
    if (!param->evt_modifiers_.empty()) {
      std::string mod_name = utils::s("Event modifier", param->evt_modifiers_.size(), false), sep;
      for (const auto& mod : param->evt_modifiers_)
        os << std::setw(wt) << mod_name << sep << utils::boldify(mod->name()) << "\n", sep = "+ ", mod_name.clear();
      os << "\n";
    }
    if (!param->out_modules_.empty()) {
      os << utils::s("Output module", param->out_modules_.size(), false);
      for (const auto& mod : param->out_modules_)
        os << "\n\t*) "
           << io::ExportModuleFactory::get().describeParameters(mod->name(), mod->parameters()).describe(1);
    }
    if (!param->taming_functions_.empty()) {
      os << std::setw(wt) << utils::s("Taming function", param->taming_functions_.size(), false) << "\n";
      for (const auto& tf : param->taming_functions_)
        os << std::setw(wt) << "" << tf->variables().at(0) << ": " << tf->expression() << "\n";
    }
    os << "\n\n"
       << std::setfill('-') << std::setw(wb + 6) << utils::boldify(" Integration parameters ") << std::setfill(' ')
       << "\n\n"
       << std::setw(wt) << "Integration" << utils::boldify(param->par_integrator.name<std::string>("N/A")) << "\n";
    for (const auto& key : param->par_integrator.keys(false))
      os << std::setw(wt) << "" << key << ": " << param->par_integrator.getString(key) << "\n";
    const auto& kin = param->process().kinematics();
    const auto& beams = kin.incomingBeams();
    os << "\n"
       << std::setfill('_') << std::setw(wb + 3) << "_/¯ EVENTS KINEMATICS ¯\\_" << std::setfill(' ') << "\n\n"
       << std::setw(wt) << "Incoming particles" << beams.positive() << ",\n"
       << std::setw(wt) << "" << beams.negative() << "\n"
       << std::setw(wt) << "C.m. energy (GeV)" << utils::format("%g", beams.sqrtS()) << "\n"
       << std::setw(wt) << "Form factors" << beams.formFactors() << "\n";
    if (beams.mode() != mode::Kinematics::ElasticElastic && !beams.structureFunctions().empty()) {
      std::ostringstream sf_list;
      for (const auto& sf : beams.structureFunctions())
        sf_list << *sf;
      os << std::setw(wt) << "Structure functions" << sf_list.str() << "\n";
    }
    os << "\n"
       << std::setfill('-') << std::setw(wb + 6) << utils::boldify(" Incoming partons ") << std::setfill(' ') << "\n\n";
    const auto& cuts = kin.cuts();
    for (const auto& lim : cuts.initial.list())  // map(particles class, limits)
      if (lim.limits.valid())
        os << std::setw(wt) << lim.description << lim.limits << "\n";
    os << "\n"
       << std::setfill('-') << std::setw(wb + 6) << utils::boldify(" Outgoing central system ") << std::setfill(' ')
       << "\n\n";
    if (!kin.minimumFinalState().empty()) {
      os << std::setw(wt) << "Minimum final state";
      std::string sep;
      for (const auto& part : kin.minimumFinalState())
        os << sep << (PDG::Id)part, sep = ", ";
      os << "\n";
    }
    for (const auto& lim : cuts.central.list())
      if (lim.limits.valid())
        os << std::setw(wt) << lim.description << lim.limits << "\n";
    if (cuts.central_particles.size() > 0) {
      os << std::setw(wt) << utils::boldify(">>> per-particle cuts:") << "\n";
      for (const auto& part_per_lim : cuts.central_particles) {
        os << " * all single " << std::setw(wt - 3) << (PDG::Id)part_per_lim.first << "\n";
        for (const auto& lim : const_cast<cuts::Central&>(part_per_lim.second).list())
          if (lim.limits.valid())
            os << "   - " << std::setw(wt - 5) << lim.description << lim.limits << "\n";
      }
    }
    os << "\n";
    os << std::setfill('-') << std::setw(wb + 6) << utils::boldify(" Proton / remnants ") << std::setfill(' ') << "\n";
    for (const auto& lim : cuts.remnants.list())
      os << "\n" << std::setw(wt) << lim.description << lim.limits;
    return os << "\n"
              << std::setfill('_') << std::setw(wb) << ""
              << "\n";
  }

  //-----------------------------------------------------------------------------------------------

  Parameters::Generation::Generation(const ParametersList& params) : SteeredObject(params) {
    (*this)
        .add("maxgen", max_gen_)
        .add("printEvery", gen_print_every_)
        .add("targetLumi", target_lumi_)
        .add("symmetrise", symmetrise_)
        .add("numThreads", num_threads_)
        .add("numPoints", num_points_);
  }

  ParametersDescription Parameters::Generation::description() {
    auto desc = ParametersDescription();
    desc.add<int>("maxgen", 0).setDescription("Number of events to generate");
    desc.add<int>("printEvery", 10000).setDescription("Printing frequency for the events content");
    desc.add<double>("targetLumi", -1.).setDescription("Target luminosity (in pb-1) to reach for this run");
    desc.add<bool>("symmetrise", false).setDescription("Are events to be symmetrised wrt beam collinear axis");
    desc.add<int>("numThreads", 1).setDescription("Number of threads to use for event generation");
    desc.add<int>("numPoints", 100);
    return desc;
  }
}  // namespace cepgen
