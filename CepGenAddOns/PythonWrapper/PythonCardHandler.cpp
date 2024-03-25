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
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"
#include "CepGenAddOns/PythonWrapper/ConfigWriter.h"

#define Py_DEBUG

namespace cepgen {
  /// CepGen Python configuration cards reader/writer
  class PythonCardHandler final : public card::Handler {
  public:
    /// Read a standard configuration card
    explicit PythonCardHandler(const ParametersList& params) : Handler(params), env_(new python::Environment(params)) {}

    PythonCardHandler& parseFile(const std::string& file) override {
      const auto filename = python::pythonPath(file);
      env_->setProgramName(filename);
      auto cfg = python::ObjectPtr::importModule(filename) /* new */;
      if (!cfg)
        throw PY_ERROR << "Failed to import the configuration card '" << filename << "'\n"
                       << " (parsed from '" << file << "').";
      parseParameters(cfg);
      parse();
      return *this;
    }
    PythonCardHandler& parseCommands(const std::vector<std::string>& str) override {
      const std::string name = "Cards.Core";
      env_->setProgramName(name);
      auto cfg = python::ObjectPtr::defineModule(name, utils::merge(str, "\n")) /* new */;
      if (!cfg)
        throw PY_ERROR << "Failed to parse a configuration string:\n"
                       << std::string(80, '-') << "\n"
                       << str << "\n"
                       << std::string(80, '-');
      parseParameters(cfg);
      parse();
      return *this;
    }

    static ParametersDescription description() {
      auto desc = Handler::description();
      desc.setDescription("Python 2/3 cards parser");
      desc.add<int>("debugging", 0).setDescription("debugging level");
      desc.add<int>("verbosity", 0).setDescription("verbosity level");
      return desc;
    }

  private:
    /// Convert the imported module into a CepGen user-steered configuration parameters object
    void parseParameters(const python::ObjectPtr& cfg) {
      CG_ASSERT(cfg);
      for (const auto& attr : cfg.attribute("__dir__")().vector<std::string>()) {
        if (attr[0] == '_')
          continue;
        const auto obj = cfg.attribute(attr);
        if (obj.is<ParametersList>())
          plist_.set(attr, obj.value<ParametersList>());
        if (obj.isVector<ParametersList>())
          plist_.set(attr, obj.vector<ParametersList>());
      }
    }
    void parse() {
      // logging module
      auto logging = plist_.get<ParametersList>("logger");
      utils::Logger::get().setLevel(logging.getAs<int, utils::Logger::Level>("level", utils::Logger::get().level()));
      utils::Logger::get().setExtended(logging.get<bool>("extended", utils::Logger::get().extended()));
      for (const auto& log_mod : logging.get<std::vector<std::string> >("enabledModules"))
        utils::Logger::get().addExceptionRule(log_mod);

      // external libraries
      for (const auto& lib : plist_.get<std::vector<std::string> >("addons"))  // additional libraries to load
        loadLibrary(lib);

      CG_DEBUG("PythonCardHandler").log([](auto& log) {
        log << "Initialised the Python cards parser.";
        for (const auto& ln : python::info())
          log << "\n\t" << ln;
      });

      // timekeeper definition (currently, does not parse the object, just check its presence)
      if (!plist_.get<ParametersList>("timer").empty())
        runParameters()->setTimeKeeper(new utils::TimeKeeper);

      // general particles definition
      if (const auto mcd_file = plist_.get<std::string>("mcdFile"); !mcd_file.empty())
        pdg::MCDFileParser::parse(mcd_file);

      // additional particles definition
      const auto parts = plist_.get<ParametersList>("PDG");
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
      if (auto process = plist_.get<ParametersList>("process"); !process.empty()) {
        process += process.get<ParametersList>("processParameters");
        process.erase("processParameters");
        auto& pkin = process.operator[]<ParametersList>("kinematics");
        {
          pkin += process.get<ParametersList>("inKinematics");
          process.erase("inKinematics");
          pkin += process.get<ParametersList>("outKinematics");
          process.erase("outKinematics");
        }
        if (process.has<int>("mode"))
          pkin.set("mode", (int)process.getAs<int, mode::Kinematics>("mode"));
        {
          if (auto& pkgen = process.operator[]<ParametersList>("kinematicsGenerator");
              pkgen.name<std::string>().empty())
            pkgen.setName<std::string>(process.get<bool>("ktFactorised", true) ? "kt2to4" : "coll2to4");
        }
        runParameters()->setProcess(ProcessFactory::get().build(process));

        for (const auto& tf : process.get<std::vector<ParametersList> >("tamingFunctions"))
          runParameters()->addTamingFunction(FunctionalFactory::get().build("ROOT", tf));
      }

      // generation parameters
      runParameters()->integrator() += plist_.get<ParametersList>("integrator");
      if (auto pgen = plist_.get<ParametersList>("generator"); !pgen.empty()) {
        runParameters()->generation().setParameters(plist_.get<ParametersList>("generator"));
        if (auto maxgen = pgen.get<int>("numEvents", -1); maxgen > 0)
          runParameters()->generation().setMaxGen(maxgen);
      }

      // event modification algorithms / hadronisers
      auto parse_evtmod_module = [&](const ParametersList& mod) {
        runParameters()->addModifier(EventModifierFactory::get().build(mod));
        auto h = runParameters()->eventModifiersSequence().rbegin()->get();
        // split the configuration into a pre-initialisation and a post-initialisation of the module parts
        h->readStrings(mod.get<std::vector<std::string> >("preConfiguration"));
        h->initialise(*runParameters());
        for (const auto& block : mod.get<std::vector<std::string> >("processConfiguration"))
          h->readStrings(mod.get<std::vector<std::string> >(block));
      };
      if (const auto had = plist_.get<ParametersList>("hadroniser"); !had.empty())  // hadronisation algorithms (legacy)
        parse_evtmod_module(had);
      for (const auto& mod :
           plist_.get<std::vector<ParametersList> >("eventSequence"))  // event modification algorithms
        parse_evtmod_module(mod);

      // output modules
      for (const auto& mod : plist_.get<std::vector<ParametersList> >("output"))
        runParameters()->addEventExporter(EventExporterFactory::get().build(mod));
    }
    void write(const std::string& filename) const override {
      python::ConfigWriter writer(ParametersList().set("filename", filename));
      writer << *runParameters();
    }
    std::unique_ptr<python::Environment> env_;
    ParametersList plist_;
  };
}  // namespace cepgen
using cepgen::PythonCardHandler;
REGISTER_CARD_HANDLER(".py", PythonCardHandler);
