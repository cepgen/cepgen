/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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

// clang-format off
#include "CepGenAddOns/PythonWrapper/Environment.h"
#include "CepGenAddOns/PythonWrapper/Error.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"
// clang-format on

#include <algorithm>

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"  // for library loading
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/TimeKeeper.h"

#define Py_DEBUG

namespace cepgen {
  /// CepGen Python configuration cards reader/writer
  class PythonCardHandler final : public card::Handler {
  public:
    /// Read a standard configuration card
    explicit PythonCardHandler(const ParametersList& params) : Handler(params), env_(new python::Environment(params)) {}

    RunParameters* parseFile(const std::string& file, RunParameters* params) override {
      const auto filename = python::pythonPath(file);
      env_->setProgramName(filename);
      if (cfg_ = python::importModule(filename) /* new */; !cfg_)
        throw PY_ERROR << "Failed to import the configuration card '" << filename << "'\n"
                       << " (parsed from '" << file << "').";
      return parse(params);
    }
    RunParameters* parseString(const std::string& str, RunParameters* params) override {
      const std::string name = "Cards.Core";
      env_->setProgramName(name);
      if (cfg_ = python::defineModule(name, str) /* new */; !cfg_)
        throw PY_ERROR << "Failed to parse a configuration string:\n"
                       << std::string(80, '-') << "\n"
                       << str << "\n"
                       << std::string(80, '-');
      return parse(params);
    }

    static ParametersDescription description() {
      auto desc = Handler::description();
      desc.setDescription("Python 2/3 cards parser");
      desc.add<int>("debugging", 0).setDescription("debugging level");
      desc.add<int>("verbosity", 0).setDescription("verbosity level");
      return desc;
    }

  private:
    RunParameters* parse(RunParameters* params) {
      CG_ASSERT(cfg_);
      // convert the imported module into a CepGen user-steered configuration parameters object
      ParametersList plist;
      for (const auto& attr : python::getVector<std::string>(python::call(cfg_.attribute("__dir__")))) {
        if (attr[0] == '_')
          continue;
        const auto obj = cfg_.attribute(attr);
        if (python::is<ParametersList>(obj))
          plist.set(attr, python::get<ParametersList>(obj));
        if (python::isVector<ParametersList>(obj))
          plist.set(attr, python::getVector<ParametersList>(obj));
      }
      rt_params_ = params;

      CG_DEBUG("PythonCardHandler").log([](auto& log) {
        log << "Initialised the Python cards parser.";
        for (const auto& ln : python::info())
          log << "\n\t" << ln;
      });

      for (const auto& lib : plist.get<std::vector<std::string> >("addons"))  // additional libraries to load
        loadLibrary(lib);

      //--- timekeeper definition
      // (currently, does not parse the object, just check its presence)
      if (!plist.get<ParametersList>("timer").empty())
        rt_params_->setTimeKeeper(new utils::TimeKeeper);

      //--- general particles definition
      if (const auto mcd_file = plist.get<std::string>("mcdFile"); !mcd_file.empty())
        pdg::MCDFileParser::parse(mcd_file);

      //--- additional particles definition
      const auto parts = plist.get<ParametersList>("PDG");
      for (const auto& k : parts.keys(true)) {
        auto props = parts.get<ParametersList>(k);
        if (props.has<int>("pdgid"))
          props.set<pdgid_t>("pdgid", props.get<int>("pdgid"));
        const ParticleProperties part(props);
        if (part.mass <= 0. && part.width <= 0.)  // skip aliases
          continue;
        if (!PDG::get().has(part.pdgid) || PDG::get()(part.pdgid) != part) {
          CG_INFO("PythonCardHandler:particles") << "Adding a new particle with PDG id=" << part.pdgid << " and name \""
                                                 << part.name << "\" to the PDG dictionary.";
          PDG::get().define(part);
        }
      }

      // process definition
      auto process = plist.get<ParametersList>("process");
      process += process.get<ParametersList>("processParameters");
      process.erase("processParameters");
      auto& pkin = process.operator[]<ParametersList>("kinematics");
      pkin += process.get<ParametersList>("inKinematics");
      process.erase("inKinematics");
      pkin += process.get<ParametersList>("outKinematics");
      process.erase("outKinematics");
      if (process.has<int>("mode"))
        pkin.set("mode", (int)process.getAs<int, mode::Kinematics>("mode"));
      rt_params_->setProcess(ProcessFactory::get().build(process));

      for (const auto& tf : process.get<std::vector<ParametersList> >("tamingFunctions"))
        rt_params_->addTamingFunction(FunctionalFactory::get().build("ROOT", tf));

      // logging module
      auto logging = plist.get<ParametersList>("logger");
      utils::Logger::get().setLevel(logging.getAs<int, utils::Logger::Level>("level", utils::Logger::get().level()));
      utils::Logger::get().setExtended(logging.get<bool>("extended", utils::Logger::get().extended()));
      for (const auto& log_mod : logging.get<std::vector<std::string> >("enabledModules"))
        utils::Logger::get().addExceptionRule(log_mod);

      auto parse_evtmod_module = [&](const ParametersList& mod) {
        rt_params_->addModifier(EventModifierFactory::get().build(mod));
        auto h = rt_params_->eventModifiersSequence().rbegin()->get();
        // split the configuration into a pre-initialisation and a post-initialisation of the module parts
        h->readStrings(mod.get<std::vector<std::string> >("preConfiguration"));
        h->initialise(*rt_params_);
        for (const auto& block : mod.get<std::vector<std::string> >("processConfiguration"))
          h->readStrings(mod.get<std::vector<std::string> >(block));
      };
      if (const auto had = plist.get<ParametersList>("hadroniser"); !had.empty())  // hadronisation algorithms (legacy)
        parse_evtmod_module(had);
      for (const auto& mod : plist.get<std::vector<ParametersList> >("eventSequence"))  // event modification algorithms
        parse_evtmod_module(mod);

      // generation parameters
      rt_params_->par_integrator += plist.get<ParametersList>("integrator");
      auto pgen = plist.get<ParametersList>("generator");
      rt_params_->generation().setParameters(pgen.set("maxgen", pgen.get<int>("numEvents")));

      for (const auto& mod : plist.get<std::vector<ParametersList> >("output"))
        rt_params_->addEventExporter(EventExporterFactory::get().build(mod));

      return rt_params_;
    }
    std::unique_ptr<python::Environment> env_;
    python::ObjectPtr cfg_{nullptr};
  };
}  // namespace cepgen
using cepgen::PythonCardHandler;
REGISTER_CARD_HANDLER(".py", PythonCardHandler);
