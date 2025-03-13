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

#include <iomanip>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/GeneratorWorkerFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

using namespace cepgen;

RunParameters::RunParameters()
    : SteeredObject(ParametersList()),
      integrator_(steer<ParametersList>("integrator")),
      generation_(steer<ParametersList>("generation")) {}

RunParameters::RunParameters(RunParameters& param)
    : SteeredObject(param),
      process_(std::move(param.process_)),
      evt_modifiers_(std::move(param.evt_modifiers_)),
      evt_exporters_(std::move(param.evt_exporters_)),
      taming_functions_(std::move(param.taming_functions_)),
      total_gen_time_(param.total_gen_time_),
      num_gen_events_(param.num_gen_events_),
      integrator_(param.integrator_),
      generation_(param.generation_),
      tmr_(std::move(param.tmr_)) {}

RunParameters::RunParameters(const RunParameters& param)
    : SteeredObject(param),
      total_gen_time_(param.total_gen_time_),
      num_gen_events_(param.num_gen_events_),
      integrator_(param.integrator_),
      generation_(param.generation_) {}

RunParameters::~RunParameters() {}  // required for unique_ptr initialisation!

RunParameters& RunParameters::operator=(RunParameters param) {
  process_ = std::move(param.process_);
  evt_modifiers_ = std::move(param.evt_modifiers_);
  evt_exporters_ = std::move(param.evt_exporters_);
  taming_functions_ = std::move(param.taming_functions_);
  total_gen_time_ = param.total_gen_time_;
  num_gen_events_ = param.num_gen_events_;
  integrator_ = param.integrator_;
  generation_ = param.generation_;
  tmr_ = std::move(param.tmr_);
  return *this;
}

void RunParameters::initialiseModules() {
  // prepare the event modifications algorithms for event generation
  for (auto& mod : evt_modifiers_)
    mod->initialise(*this);
  // prepare the output modules for event generation
  for (auto& mod : evt_exporters_)
    mod->initialise(*this);
}

void RunParameters::prepareRun() {
  if (tmr_)
    tmr_->clear();
  CG_TICKER(tmr_.get());

  //--- clear the run statistics
  total_gen_time_ = 0.;
  num_gen_events_ = 0ul;
}

void RunParameters::setTimeKeeper(utils::TimeKeeper* kpr) { tmr_.reset(kpr); }

void RunParameters::addGenerationTime(double gen_time) {
  total_gen_time_ += gen_time;
  num_gen_events_++;
}

proc::Process& RunParameters::process() {
  if (!process_)
    throw CG_FATAL("RunParameters:process") << "Failed to retrieve a process configuration block.";
  return *process_.get();
}

const proc::Process& RunParameters::process() const {
  if (!process_)
    throw CG_FATAL("RunParameters:process") << "Failed to retrieve a process configuration block.";
  return *process_.get();
}

std::string RunParameters::processName() const {
  if (!process_)
    return "no process";
  return process_->name();
}

void RunParameters::clearProcess() { delete process_.release(); }

void RunParameters::setProcess(std::unique_ptr<proc::Process> proc) { process_ = std::move(proc); }

void RunParameters::setProcess(proc::Process* proc) {
  if (!proc)
    throw CG_FATAL("RunParameters") << "Trying to clone an invalid process!";
  process_.reset(proc);
}

const Kinematics& RunParameters::kinematics() const {
  if (!process_)
    throw CG_FATAL("RunParameters") << "Process must be defined before its kinematics is retrieved!";
  return process_->kinematics();
}

EventModifier& RunParameters::eventModifier(size_t i) const { return *evt_modifiers_.at(i); }

void RunParameters::clearEventModifiersSequence() { evt_modifiers_.clear(); }

void RunParameters::addModifier(std::unique_ptr<EventModifier> mod) { evt_modifiers_.emplace_back(std::move(mod)); }

void RunParameters::addModifier(EventModifier* mod) {
  evt_modifiers_.emplace_back(std::unique_ptr<EventModifier>(mod));
}

EventExporter& RunParameters::eventExporter(size_t i) const { return *evt_exporters_.at(i); }

void RunParameters::clearEventExportersSequence() { evt_exporters_.clear(); }

void RunParameters::addEventExporter(std::unique_ptr<EventExporter> mod) {
  evt_exporters_.emplace_back(std::move(mod));
}

void RunParameters::addEventExporter(EventExporter* mod) {
  evt_exporters_.emplace_back(std::unique_ptr<EventExporter>(mod));
}

void RunParameters::addTamingFunction(std::unique_ptr<utils::Functional> fct) {
  taming_functions_.emplace_back(std::move(fct));
}

ParametersDescription RunParameters::description() {
  auto desc = ParametersDescription();
  desc.add<ParametersDescription>("integrator", IntegratorFactory::get().describeParameters("Vegas"));
  desc.add<ParametersDescription>("generation", Generation::description());
  return desc;
}

namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const RunParameters& param) {
    constexpr int wb = 90, wt = 33;

    os << std::left << "\n"
       << std::setfill('_') << std::setw(wb + 3) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill(' ') << "\n\n";
    if (param.process_) {
      const auto& proc_params = param.process().parameters();
      os << std::setw(wt) << "Process to generate:"
         << utils::boldify(ProcessFactory::get().describeParameters(proc_params).description()) << "\n";
      for (const auto& key : proc_params.keys(false)) {
        if (key == "kinematics" || key == "partonFluxes" || key == "ktFluxes")  // these are shown below
          continue;
        os << std::setw(wt) << "" << key << ": "
           << (proc_params.has<ParametersList>(key) ? proc_params.get<ParametersList>(key).print(true)
                                                    : proc_params.getString(key))
           << "\n";
      }
    }
    if (!param.evt_modifiers_.empty() || !param.evt_exporters_.empty() || !param.taming_functions_.empty())
      os << "\n"
         << std::setfill('-') << std::setw(wb + 6) << utils::boldify(" Event treatment ") << std::setfill(' ')
         << "\n\n";
    if (!param.evt_modifiers_.empty()) {
      std::string mod_name = utils::s("Event modifier", param.evt_modifiers_.size(), false), sep;
      for (const auto& mod : param.evt_modifiers_)
        os << std::setw(wt) << mod_name << sep << utils::boldify(mod->name()) << "\n", sep = "+ ", mod_name.clear();
      os << "\n";
    }
    if (!param.evt_exporters_.empty()) {
      os << utils::s("Output module", param.evt_exporters_.size(), false);
      for (const auto& mod : param.evt_exporters_)
        os << "\n\t*) " << EventExporterFactory::get().describeParameters(mod->name(), mod->parameters()).describe(1);
    }
    if (!param.taming_functions_.empty()) {
      os << std::setw(wt) << utils::s("Taming function", param.taming_functions_.size(), false) << "\n";
      for (const auto& tf : param.taming_functions_)
        os << std::setw(wt) << "" << tf->variables().at(0) << ": " << tf->expression() << "\n";
    }
    os << "\n\n"
       << std::setfill('-') << std::setw(wb + 6) << utils::boldify(" Integration/generation parameters ")
       << std::setfill(' ') << "\n\n"
       << std::setw(wt) << "Integration" << utils::boldify(param.integrator_.name("N/A")) << "\n";
    for (const auto& key : param.integrator_.keys(false))
      os << std::setw(wt) << "" << key << ": " << param.integrator_.getString(key) << "\n";
    os << std::setw(wt) << "Event generation? " << utils::yesno(param.generation_.enabled()) << "\n"
       << std::setw(wt) << "Number of events to generate" << utils::boldify(param.generation_.maxGen()) << "\n"
       << std::setw(wt) << "Generator worker"
       << param.generation_.parameters().get<ParametersList>("worker").print(true) << "\n";
    if (param.generation_.numThreads() > 1)
      os << std::setw(wt) << "Number of threads" << param.generation_.numThreads() << "\n";
    os << std::setw(wt) << "Number of points to try per bin" << param.generation_.numPoints() << "\n"
       << std::setw(wt) << "Verbosity level " << utils::Logger::get().level() << "\n";
    const auto& kin = param.process().kinematics();
    const auto& beams = kin.incomingBeams();
    os << "\n"
       << std::setfill('_') << std::setw(wb + 3) << "_/¯ EVENTS KINEMATICS ¯\\_" << std::setfill(' ') << "\n\n"
       << std::setw(wt) << "Incoming particles" << beams.positive() << ",\n"
       << std::setw(wt) << "" << beams.negative() << "\n"
       << std::setw(wt) << "C.m. energy (GeV)" << utils::format("%g", beams.sqrtS()) << "\n";
    if (beams.mode() != mode::Kinematics::ElasticElastic)
      os << std::setw(wt) << "Structure functions"
         << utils::boldify(
                StructureFunctionsFactory::get().describeParameters(beams.structureFunctions()).description())
         << "\n"
         << std::setw(wt) << "" << beams.structureFunctions().print(true) << "\n ";

    // helper to print cuts list for a given category
    const auto dump_cuts = [&os](const auto& obj) {
      for (const auto& lim : obj.parameters().template keysOf<Limits>())
        if (const auto& limit = obj.parameters().template get<Limits>(lim); limit.valid() && obj.description().has(lim))
          os << std::setw(wt) << obj.description().get(lim).description() << limit << "\n";
      for (const auto& vlim : obj.parameters().template keysOf<std::vector<Limits> >())
        if (const auto& limit = obj.parameters().template get<std::vector<Limits> >(vlim); obj.description().has(vlim))
          os << std::setw(wt) << obj.description().get(vlim).description() << utils::repr(limit, " and ") << "\n";
    };

    os << "\n"
       << std::setfill('-') << std::setw(wb + 6) << utils::boldify(" Incoming partons ") << std::setfill(' ') << "\n\n";
    const auto& cuts = kin.cuts();
    dump_cuts(cuts.initial);
    os << "\n"
       << std::setfill('-') << std::setw(wb + 6) << utils::boldify(" Outgoing central system ") << std::setfill(' ')
       << "\n\n";
    if (!kin.minimumFinalState().empty()) {
      os << std::setw(wt) << "Minimum final state";
      std::string sep;
      for (const auto& part : kin.minimumFinalState())
        os << sep << static_cast<PDG::Id>(part), sep = ", ";
      os << "\n";
    }
    dump_cuts(cuts.central);
    if (cuts.central_particles.size() > 0) {
      os << std::setw(wt) << utils::boldify(">>> per-particle cuts:") << "\n";
      for (const auto& part_per_lim : cuts.central_particles) {
        os << " * all single " << std::setw(wt - 3) << static_cast<PDG::Id>(part_per_lim.first) << "\n";
        for (const auto& lim : part_per_lim.second.parameters().keysOf<Limits>()) {
          const auto& limit = part_per_lim.second.parameters().get<Limits>(lim);
          if (limit.valid())
            os << "   - " << std::setw(wt - 5) << cuts::Central::description().get(lim).description() << limit << "\n";
        }
      }
    }
    os << "\n";
    os << std::setfill('-') << std::setw(wb + 6) << utils::boldify(" Proton / remnants ") << std::setfill(' ')
       << "\n\n";
    dump_cuts(cuts.remnants);
    return os << "\n"
              << std::setfill('_') << std::setw(wb) << ""
              << "\n";
  }
}  // namespace cepgen

//-----------------------------------------------------------------------------------------------

RunParameters::Generation::Generation(const ParametersList& params) : SteeredObject(params) {
  (*this)
      .add("maxgen", max_gen_)
      .add("printEvery", gen_print_every_)
      .add("targetLumi", target_lumi_)
      .add("symmetrise", symmetrise_)
      .add("numThreads", num_threads_)
      .add("numPoints", num_points_);
}

ParametersDescription RunParameters::Generation::description() {
  auto desc = ParametersDescription();
  desc.add<ParametersDescription>("worker", GeneratorWorkerFactory::get().describeParameters("grid_optimised"))
      .setDescription("type of generator worker to use for event generation");
  desc.add<int>("maxgen", 0).setDescription("Number of events to generate");
  desc.add<int>("printEvery", 10000).setDescription("Printing frequency for the events content");
  desc.add<double>("targetLumi", -1.).setDescription("Target luminosity (in pb-1) to reach for this run");
  desc.add<bool>("symmetrise", false).setDescription("Are events to be symmetrised wrt beam collinear axis");
  desc.add<int>("numThreads", 1).setDescription("Number of threads to use for event generation");
  desc.add<int>("numPoints", 100);
  return desc;
}
