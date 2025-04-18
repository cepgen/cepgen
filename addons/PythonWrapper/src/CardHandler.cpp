/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2025  Laurent Forthomme
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
#include "CepGenPython/Environment.h"
#include "CepGenPython/Error.h"
#include "CepGenPython/Utils.h"
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
#include "CepGenPython/ConfigWriter.h"

#define Py_DEBUG

using namespace std::string_literals;

namespace cepgen::python {
  /// CepGen Python configuration cards reader/writer
  class CardHandler final : public card::Handler {
  public:
    /// Read a standard configuration card
    explicit CardHandler(const ParametersList& params)
        : Handler(params), plist_(params_.operator[]<ParametersList>("parsed")) {}

    CardHandler& parseFile(const std::string& file) override {
      const auto filename = pythonPath(file);
      Environment environment(params_);
      environment.setProgramName(filename);
      if (const auto cfg = ObjectPtr::importModule(filename) /* new */; cfg)
        parseParameters(cfg);
      else
        throw PY_ERROR << "Failed to import the configuration card '" << filename << "'\n"
                       << " (parsed from '" << file << "').";
      parse();
      return *this;
    }
    CardHandler& parseCommands(const std::vector<std::string>& str) override {
      const std::string name = "Cards.Core";
      Environment environment(params_);
      environment.setProgramName(name);
      if (const auto cfg = ObjectPtr::defineModule(name, utils::merge(str, "\n")) /* new */; cfg)
        parseParameters(cfg);
      else
        throw PY_ERROR << "Failed to parse a configuration string:\n"
                       << std::string(80, '-') << "\n"
                       << str << "\n"
                       << std::string(80, '-');
      parse();
      return *this;
    }

    static ParametersDescription description() {
      auto desc = Handler::description();
      desc.setDescription("Python 2/3 cards parser");
      desc.add("debugging", 0).setDescription("debugging level");
      desc.add("verbosity", 0).setDescription("verbosity level");
      return desc;
    }

  private:
    /// Convert the imported module into a CepGen user-steered configuration parameters object
    void parseParameters(const ObjectPtr& cfg) const {
      CG_ASSERT(cfg);
      for (const auto& attr : cfg.attribute("__dir__")().vector<std::string>()) {
        if (attr[0] == '_')
          continue;
        const auto obj = cfg.attribute(attr);
        if (obj.is<ParametersList>())
          plist_.set(attr, obj.value<ParametersList>());
        else if (obj.isVector<ParametersList>())
          plist_.set(attr, obj.vector<ParametersList>());
        else if (obj.is<int>())
          plist_.set(attr, obj.value<int>());
        else if (obj.isVector<int>())
          plist_.set(attr, obj.vector<int>());
        else if (obj.is<double>())
          plist_.set(attr, obj.value<double>());
        else if (obj.isVector<double>()) {
          if (obj.is<Limits>())
            plist_.set(attr, obj.value<Limits>());
          plist_.set(attr, obj.vector<double>());
        } else if (obj.is<std::string>())
          plist_.set(attr, obj.value<std::string>());
        else if (obj.isVector<std::string>())
          plist_.set(attr, obj.vector<std::string>());
        else if (obj.isVector<Limits>())
          plist_.set(attr, obj.vector<Limits>());
      }
    }
    void parse() {
      // logging module
      const auto logging = plist_.get<ParametersList>("logger");
      utils::Logger::get().setLevel(logging.getAs<int, utils::Logger::Level>("level", utils::Logger::get().level()));
      utils::Logger::get().setExtended(logging.get<bool>("extended", utils::Logger::get().extended()));
      for (const auto& log_mod : logging.get<std::vector<std::string> >("enabledModules"))
        utils::Logger::get().addExceptionRule(log_mod);

      // external libraries
      for (const auto& lib : plist_.get<std::vector<std::string> >("addons"))  // additional libraries to load
        loadLibrary(lib);

      CG_DEBUG("python:CardHandler").log([](auto& log) {
        log << "Initialised the Python cards parser.";
        for (const auto& ln : info())
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
          CG_INFO("python:CardHandler:particles") << "Adding a new particle with PDG id=" << part.pdgid
                                                  << " and name \"" << part.name << "\" to the PDG dictionary.";
          PDG::get().define(part);
        }
      }

      // process definition
      if (auto process = plist_.get<ParametersList>("process"); !process.empty()) {
        auto& kinematics_parameters = process.operator[]<ParametersList>("kinematics");
        {  // remove extra layer of 'processParameters' and move it to the main process parameters block
          process += process.get<ParametersList>("processParameters");
          process.erase("processParameters");
          if (process.has<int>("mode")) {  // move the kinematics mode from process to the main kinematics block
            kinematics_parameters.set("mode", static_cast<int>(process.getAs<int, mode::Kinematics>("mode")));
            process.erase("mode");
          }
        }
        {  // remove extra layers of 'inKinematics' and 'outKinematics' and move them to the main kinematics block
          kinematics_parameters += process.get<ParametersList>("inKinematics");
          process.erase("inKinematics");
          kinematics_parameters += process.get<ParametersList>("outKinematics");
          process.erase("outKinematics");
        }
        if (auto& pkgen = process.operator[]<ParametersList>("kinematicsGenerator"); pkgen.name().empty())
          pkgen.setName((process.get<bool>("ktFactorised", true) ? "kt"s : "coll"s) + ":2to4");
        runParameters()->setProcess(ProcessFactory::get().build(process));

        for (const auto& tf : process.get<std::vector<ParametersList> >("tamingFunctions"))
          runParameters()->addTamingFunction(FunctionalFactory::get().build("python", tf));
      }

      // generation parameters
      runParameters()->integrator() += plist_.get<ParametersList>("integrator");
      if (const auto pgen = plist_.get<ParametersList>("generator"); !pgen.empty()) {
        runParameters()->generation().setParameters(plist_.get<ParametersList>("generator"));
        if (const auto num_events = pgen.get<int>("numEvents", -1); num_events > 0)
          runParameters()->generation().setMaxGen(num_events);
      }

      // event modification algorithms / hadronisers
      auto parse_event_modifier_parameters = [&](const ParametersList& mod) {
        runParameters()->addModifier(EventModifierFactory::get().build(mod));
        const auto h = runParameters()->eventModifiersSequence().rbegin()->get();
        // split the configuration into a pre-initialisation and a post-initialisation of the module parts
        h->readStrings(mod.get<std::vector<std::string> >("preConfiguration"));
        h->initialise(*runParameters());
        for (const auto& block : mod.get<std::vector<std::string> >("processConfiguration"))
          h->readStrings(mod.get<std::vector<std::string> >(block));
      };
      if (const auto had = plist_.get<ParametersList>("hadroniser"); !had.empty())  // hadronisation algorithms (legacy)
        parse_event_modifier_parameters(had);
      for (const auto& mod :
           plist_.get<std::vector<ParametersList> >("eventSequence"))  // event modification algorithms
        parse_event_modifier_parameters(mod);

      // output modules
      for (const auto& mod : plist_.get<std::vector<ParametersList> >("output"))
        runParameters()->addEventExporter(EventExporterFactory::get().build(mod));
    }
    void write(const std::string& filename) const override {
      ConfigWriter writer(ParametersList().set("filename", filename));
      writer << *runParameters();
    }
    ParametersList& plist_;
  };
}  // namespace cepgen::python
using PythonCardHandler = cepgen::python::CardHandler;
REGISTER_CARD_HANDLER(".py", PythonCardHandler);
