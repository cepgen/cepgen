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

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

using namespace cepgen;
using namespace cepgen::card;
using namespace std::string_literals;

/// Command line parser
class CommandLineHandler final : public Handler {
public:
  /// Cast command line arguments into a configuration word
  explicit CommandLineHandler(const ParametersList& params) : Handler(params) {
    if (const auto argv = steer<std::vector<std::string> >("args"); !argv.empty())
      parseCommands(argv);
  }

  static ParametersDescription description() {
    auto desc = Handler::description();
    desc.setDescription("Command line configuration parser");
    desc.add("args", std::vector<std::string>{}).setDescription("Collection of arguments to be parsed");
    return desc;
  }

  CommandLineHandler& parseFile(const std::string& filename) override {
    if (filename.empty())
      throw CG_FATAL("CommandLineHandler") << "Empty filename to be parsed! Aborting.";
    return parseCommands(utils::split(utils::readFile(filename), '\n'));
  }

  CommandLineHandler& parseCommands(const std::vector<std::string>& commands) override {
    if (commands.empty())
      return *this;
    ParametersList parameters;
    for (const auto& command : commands)
      parameters.feed(command);
    CG_INFO("CommandLineHandler") << "Arguments list: " << commands << " unpacked to:\n\t" << parameters << ".";

    // timer definition
    if (parameters.get<bool>("timer"s, false))
      runParameters()->setTimeKeeper(new utils::TimeKeeper);

    // logging definition
    if (parameters.has<int>("logging"s))
      utils::Logger::get().setLevel(parameters.getAs<int, utils::Logger::Level>("logging"s));
    else if (parameters.has<ParametersList>("logging"s)) {
      const auto& logging_params = parameters.get<ParametersList>("logging"s);
      if (logging_params.has<int>("level"s))
        utils::Logger::get().setLevel(logging_params.getAs<int, utils::Logger::Level>("level"s));
      if (logging_params.has<std::string>("modules"s))
        utils::Logger::get().addExceptionRule(logging_params.get<std::string>("modules"s));
      else if (logging_params.has<std::vector<std::string> >("modules"s))
        for (const auto& mod : logging_params.get<std::vector<std::string> >("modules"s))
          utils::Logger::get().addExceptionRule(mod);
      utils::Logger::get().setExtended(logging_params.get<bool>("extended"s, false));
    }

    // PDG definition
    const auto pdg_params = parameters.get<ParametersList>("pdg"s);
    for (const auto& id : pdg_params.keys())
      PDG::get().define(pdg_params.get<ParticleProperties>(id));

    // phase space definition
    auto kinematics_params = parameters.get<ParametersList>("kinematics"s);

    // process definition
    if (auto process_params = parameters.get<ParametersList>("process"s); !process_params.empty()) {
      if (runParameters()->hasProcess())
        process_params = ParametersList(runParameters()->process().parameters()) + process_params;
      runParameters()->setProcess(ProcessFactory::get().build(process_params));
      if (process_params.has<int>("mode"s))
        kinematics_params.set("mode"s, process_params.get<int>("mode"s));
    }

    if (!kinematics_params.empty()) {  // set auxiliary information for phase space definition
      if (kinematics_params.has<int>("strfun"))
        kinematics_params
            .set(
                "structureFunctions"s,
                StructureFunctionsFactory::get().describeParameters(kinematics_params.get<int>("strfun"s)).parameters())
            .erase("strfun"s);
      else if (kinematics_params.has<std::string>("strfun"s))
        kinematics_params
            .set("structureFunctions"s,
                 StructureFunctionsFactory::get()
                     .describeParameters(kinematics_params.get<std::string>("strfun"s))
                     .parameters())
            .erase("strfun"s);
      else if (kinematics_params.has<ParametersList>("strfun"s))
        kinematics_params.rename("strfun"s, "structureFunctions");
      kinematics_params.rename("formfac"s, "formFactors");

      // get the kinematics as already defined in the process object and modify it accordingly
      kinematics_params = runParameters()->process().kinematics().parameters() + kinematics_params;
      runParameters()->process().kinematics().setParameters(kinematics_params);
    }

    // integration
    parameters.fill<ParametersList>("integrator"s, runParameters()->integrator());

    // events generation
    const auto& generation_params = parameters.get<ParametersList>("generation"s);
    runParameters()->generation().setMaxGen(
        generation_params.get<int>("ngen"s, runParameters()->generation().maxGen()));
    if (generation_params.has<int>("nthreads"s))
      runParameters()->generation().setNumThreads(generation_params.get<int>("nthreads"s));
    if (generation_params.has<int>("nprn"s))
      runParameters()->generation().setPrintEvery(generation_params.get<int>("nprn"s));
    if (generation_params.has<int>("seed"s))
      runParameters()->integrator().set("seed", generation_params.get<int>("seed"s));

    // event modification modules
    if (const auto& event_mod_params = parameters.get<ParametersList>("eventmod"s);
        !event_mod_params.keys(true).empty()) {
      runParameters()->addModifier(EventModifierFactory::get().build(event_mod_params));
      runParameters()->eventModifiersSequence().rbegin()->get()->initialise(*runParameters());
    }

    // output modules definition
    if (const auto& output_params = parameters.get<ParametersList>("output"s); !output_params.keys(true).empty())
      runParameters()->addEventExporter(EventExporterFactory::get().build(output_params));
    return *this;
  }
};
REGISTER_CARD_HANDLER(".cmd", CommandLineHandler);
