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
    desc.add<std::vector<std::string> >("args", {}).setDescription("Collection of arguments to be parsed");
    return desc;
  }

  inline CommandLineHandler& parseFile(const std::string& filename) override {
    if (filename.empty())
      throw CG_FATAL("CommandLineHandler") << "Empty filename to be parsed! Aborting.";
    return parseCommands(utils::split(utils::readFile(filename), '\n'));
  }

  inline CommandLineHandler& parseCommands(const std::vector<std::string>& commands) override {
    if (commands.empty())
      return *this;
    ParametersList pars;
    for (const auto& cmd : commands)
      pars.feed(cmd);
    CG_INFO("CommandLineHandler") << "Arguments list: " << commands << " unpacked to:\n\t" << pars << ".";

    //----- timer definition
    if (pars.get<bool>("timer", false))
      runParameters()->setTimeKeeper(new utils::TimeKeeper);

    //----- logging definition
    if (pars.has<int>("logging"))
      utils::Logger::get().setLevel(pars.getAs<int, cepgen::utils::Logger::Level>("logging"));
    else if (pars.has<ParametersList>("logging")) {
      const auto& log = pars.get<ParametersList>("logging");
      if (log.has<int>("level"))
        utils::Logger::get().setLevel(log.getAs<int, cepgen::utils::Logger::Level>("level"));
      if (log.has<std::string>("modules"))
        utils::Logger::get().addExceptionRule(log.get<std::string>("modules"));
      else if (log.has<std::vector<std::string> >("modules"))
        for (const auto& mod : log.get<std::vector<std::string> >("modules"))
          utils::Logger::get().addExceptionRule(mod);
      utils::Logger::get().setExtended(log.get<bool>("extended", false));
    }

    //----- PDG definition
    auto pars_pdg = pars.get<ParametersList>("pdg");
    for (const auto& id : pars_pdg.keys())
      PDG::get().define(pars_pdg.get<ParticleProperties>(id));

    //----- phase space definition
    auto pars_kin = pars.get<ParametersList>("kinematics");

    //----- process definition
    if (auto proc = pars.get<ParametersList>("process"); !proc.empty()) {
      if (runParameters()->hasProcess())
        proc = ParametersList(runParameters()->process().parameters()) + proc;
      runParameters()->setProcess(ProcessFactory::get().build(proc));
      if (proc.has<int>("mode"))
        pars_kin.set<int>("mode", proc.get<int>("mode"));
    }

    if (!pars_kin.empty()) {
      //----- set auxiliary information for phase space definition
      if (pars_kin.has<int>("strfun"))
        pars_kin
            .set<ParametersList>(
                "structureFunctions",
                StructureFunctionsFactory::get().describeParameters(pars_kin.get<int>("strfun")).parameters())
            .erase("strfun");
      else if (pars_kin.has<std::string>("strfun"))
        pars_kin
            .set<ParametersList>(
                "structureFunctions",
                StructureFunctionsFactory::get().describeParameters(pars_kin.get<std::string>("strfun")).parameters())
            .erase("strfun");
      else if (pars_kin.has<ParametersList>("strfun"))
        pars_kin.rename("strfun", "structureFunctions");
      pars_kin.rename("formfac", "formFactors");

      //----- get the kinematics as already defined in the process object and modify it accordingly
      pars_kin = runParameters()->process().kinematics().parameters() + pars_kin;
      runParameters()->process().kinematics().setParameters(pars_kin);
    }

    //----- integration
    pars.fill<ParametersList>("integrator", runParameters()->integrator());

    //----- events generation
    const auto& gen = pars.get<ParametersList>("generation");
    runParameters()->generation().setMaxGen(gen.get<int>("ngen", runParameters()->generation().maxGen()));
    if (gen.has<int>("nthreads"))
      runParameters()->generation().setNumThreads(gen.get<int>("nthreads"));
    if (gen.has<int>("nprn"))
      runParameters()->generation().setPrintEvery(gen.get<int>("nprn"));
    if (gen.has<int>("seed"))
      runParameters()->integrator().set<int>("seed", gen.get<int>("seed"));

    //----- event modification modules
    if (const auto& mod = pars.get<ParametersList>("eventmod"); !mod.keys(true).empty()) {
      runParameters()->addModifier(EventModifierFactory::get().build(mod));
      runParameters()->eventModifiersSequence().rbegin()->get()->initialise(*runParameters());
    }

    //----- output modules definition
    if (const auto& out = pars.get<ParametersList>("output"); !out.keys(true).empty())
      runParameters()->addEventExporter(EventExporterFactory::get().build(out));
    return *this;
  }
};
REGISTER_CARD_HANDLER(".cmd", CommandLineHandler);
