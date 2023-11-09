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

// clang-format off
#include "CepGenAddOns/PythonWrapper/Environment.h"
#include "CepGenAddOns/PythonWrapper/Error.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"
// clang-format on

#include <algorithm>

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"  // for library loading
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/TimeKeeper.h"

#define Py_DEBUG

namespace cepgen {
  namespace card {
    /// CepGen Python configuration cards reader/writer
    class PythonHandler final : public Handler {
    public:
      /// Read a standard configuration card
      explicit PythonHandler(const ParametersList&);
      Parameters* parseFile(const std::string&, Parameters*) override;
      Parameters* parseString(const std::string&, Parameters*) override;

      static ParametersDescription description();

    private:
      static constexpr const char* ADDONS_NAME = "addons";
      static constexpr const char* TIMER_NAME = "timer";
      static constexpr const char* PROCESS_NAME = "process";
      static constexpr const char* HADR_NAME = "hadroniser";
      static constexpr const char* EVT_MOD_SEQ_NAME = "eventSequence";
      static constexpr const char* LOGGER_NAME = "logger";
      static constexpr const char* INTEGRATOR_NAME = "integrator";
      static constexpr const char* GENERATOR_NAME = "generator";
      static constexpr const char* OUTPUT_NAME = "output";

      static constexpr const char* PDGLIST_NAME = "PDG";
      static constexpr const char* MCD_NAME = "mcdFile";

      Parameters* parse(Parameters*);
      void parseLogging(PyObject*);
      void parseIntegrator(PyObject*);
      void parseGenerator(PyObject*);
      void parseHadroniser(PyObject*);
      void parseEventModifiers(PyObject*);
      void parseEventExporter(PyObject*);
      void parseEventExporters(PyObject*);
      void parseExtraParticles(PyObject*);

      std::unique_ptr<python::Environment> env_;
      python::ObjectPtr cfg_{nullptr};
    };

    PythonHandler::PythonHandler(const ParametersList& params)
        : Handler(params), env_(new python::Environment(params)) {}

    Parameters* PythonHandler::parseFile(const std::string& file, Parameters* params) {
      std::string filename = python::pythonPath(file);
      env_->setProgramName(filename);
      cfg_ = python::importModule(filename);  // new
      if (!cfg_)
        throw PY_ERROR << "Failed to import the configuration card '" << filename << "'\n"
                       << " (parsed from '" << file << "').";
      return parse(params);
    }

    Parameters* PythonHandler::parseString(const std::string& str, Parameters* params) {
      env_->setProgramName("Cards.Core");
      cfg_ = python::defineModule("Cards.Core", str);  // new
      if (!cfg_)
        throw PY_ERROR << "Failed to parse a configuration string:\n"
                       << std::string(80, '-') << "\n"
                       << str << "\n"
                       << std::string(80, '-');
      return parse(params);
    }

    Parameters* PythonHandler::parse(Parameters* params) {
      if (!cfg_)
        throw PY_ERROR << "Python configuration card was not defined.";

      rt_params_ = params;

      CG_DEBUG("PythonHandler").log([](auto& log) {
        log << "Initialised the Python cards parser.";
        for (const auto& ln : python::info())
          log << "\n\t" << ln;
      });

      auto parseAttr = [this](const std::string& name, std::function<void(PyObject*)> callback) -> void {
        auto pobj = python::getAttribute(cfg_.get(), name);
        if (pobj)
          callback(pobj.get());
      };

      //--- additional libraries to load
      parseAttr(ADDONS_NAME, [](PyObject* padd) {
        for (const auto& lib : python::getVector<std::string>(padd))
          loadLibrary(lib);
      });

      //--- timekeeper definition
      // (currently, does not parse the object, just check its presence)
      parseAttr(TIMER_NAME, [this](PyObject*) { rt_params_->setTimeKeeper(new utils::TimeKeeper); });

      //--- general particles definition
      parseAttr(MCD_NAME, [](PyObject* ppdg) { pdg::MCDFileParser::parse(python::get<std::string>(ppdg).c_str()); });

      //--- additional particles definition
      parseAttr(PDGLIST_NAME, [this](PyObject* pextp) { parseExtraParticles(pextp); });

      //--- process definition
      parseAttr(PROCESS_NAME, [this](PyObject* process) {
        //--- list of process-specific parameters
        ParametersList proc_params;
        python::fillParameter(process, "processParameters", proc_params);

        //--- type of process to consider
        auto* pproc_name = python::element(process, MODULE_NAME);  // borrowed
        if (!pproc_name)
          PY_ERROR << "Failed to extract the process name from the configuration.";

        //--- process mode
        const auto proc_name = python::get<std::string>(pproc_name);
        if (auto* pkt = python::element(process, "ktFactorised"))
          proc_params.set<bool>("ktFactorised", python::get<bool>(pkt));
        CG_DEBUG("PythonHandler") << "Building a process with name '" << proc_name << "' and parameters:\n\t"
                                  << proc_params << ".";
        //--- process kinematics
        auto pkin = ParametersList();
        if (auto* pin_kinematics = python::element(process, "inKinematics"))  // borrowed
          pkin += python::get<ParametersList>(pin_kinematics);

        if (auto* pout_kinematics = python::element(process, "outKinematics"))  // borrowed
          pkin += python::get<ParametersList>(pout_kinematics);

        if (proc_params.has<int>("mode"))
          pkin.set<int>("mode", proc_params.get<int>("mode"));
        CG_DEBUG("PythonHandler") << "Setting kinematics to:\n" << pkin << ".";

        auto proc_obj = ProcessFactory::get().build(proc_name, proc_params.set("kinematics", pkin));

        // feed the runtime parameters with this newly populated process object
        rt_params_->setProcess(std::move(proc_obj));

        //--- taming functions
        if (auto* ptam = python::element(process, "tamingFunctions"))  // borrowed
          for (const auto& p : python::getVector<ParametersList>(ptam))
            rt_params_->addTamingFunction(FunctionalFactory::get().build("ROOT", p));
      });

      parseAttr(LOGGER_NAME, [this](PyObject* plog) { parseLogging(plog); });

      //--- hadroniser parameters (legacy)
      parseAttr(HADR_NAME, [this](PyObject* phad) { parseHadroniser(phad); });
      parseAttr(EVT_MOD_SEQ_NAME, [this](PyObject* pmod_seq) { parseEventModifiers(pmod_seq); });

      //--- generation parameters
      parseAttr(INTEGRATOR_NAME, [this](PyObject* pint) { parseIntegrator(pint); });
      parseAttr(GENERATOR_NAME, [this](PyObject* pgen) { parseGenerator(pgen); });
      parseAttr(OUTPUT_NAME, [this](PyObject* pout) { parseEventExporters(pout); });

      return rt_params_;
    }

    void PythonHandler::parseLogging(PyObject* log) {
      int log_level{(int)utils::Logger::get().level()};
      python::fillParameter(log, "level", log_level);
      utils::Logger::get().setLevel((utils::Logger::Level)log_level);
      bool extended{utils::Logger::get().extended()};
      python::fillParameter(log, "extended", extended);
      utils::Logger::get().setExtended(extended);
      std::vector<std::string> enabled_modules;
      python::fillParameter(log, "enabledModules", enabled_modules);
      for (const auto& mod : enabled_modules)
        utils::Logger::get().addExceptionRule(mod);
    }

    void PythonHandler::parseIntegrator(PyObject* integr) {
      if (!python::is<ParametersList>(integr))
        PY_ERROR << "Integrator object should be a dictionary.";
      rt_params_->par_integrator += python::get<ParametersList>(integr);
    }

    void PythonHandler::parseGenerator(PyObject* gen) {
      if (!python::is<ParametersList>(gen))
        PY_ERROR << "Generation information object should be a dictionary.";
      auto plist = python::get<ParametersList>(gen);
      plist.set<int>("maxgen", plist.get<int>("numEvents"));
      rt_params_->generation().setParameters(plist);
    }

    void PythonHandler::parseEventModifiers(PyObject* mod) {
      if (!python::isVector<ParametersList>(mod))
        PY_ERROR << "Event modification definition object should be a list/Sequence.";

      for (Py_ssize_t i = 0; i < PyList_Size(mod); ++i)
        parseHadroniser(PyList_GetItem(mod, i));
    }

    void PythonHandler::parseHadroniser(PyObject* mod) {
      if (!python::is<ParametersList>(mod))
        PY_ERROR << "Event modification definition object should be a dictionary.";

      auto* pname = python::element(mod, MODULE_NAME);  // borrowed
      if (!pname)
        PY_ERROR << "Event modification algorithm name is required.";
      std::string mod_name = python::get<std::string>(pname);

      rt_params_->addModifier(EventModifierFactory::get().build(mod_name, python::get<ParametersList>(mod)));

      auto h = rt_params_->eventModifiersSequence().rbegin()->get();
      // split the configuration into a pre-initialisation and a post-initialisation of the module parts
      {
        std::vector<std::string> config;
        python::fillParameter(mod, "preConfiguration", config);
        h->readStrings(config);
      }
      h->initialise(*rt_params_);
      {
        std::vector<std::string> config;
        python::fillParameter(mod, "processConfiguration", config);
        for (const auto& block : config) {
          std::vector<std::string> config_blk;
          python::fillParameter(mod, block.c_str(), config_blk);
          h->readStrings(config_blk);
        }
      }
    }

    void PythonHandler::parseEventExporters(PyObject* mod) {
      if (!python::isVector<ParametersList>(mod))
        PY_ERROR << "Output modules definition object should be a list/Sequence.";

      for (Py_ssize_t i = 0; i < PyList_Size(mod); ++i)
        parseEventExporter(PyList_GetItem(mod, i));
    }

    void PythonHandler::parseEventExporter(PyObject* pout) {
      if (!python::is<ParametersList>(pout))
        PY_ERROR << "Invalid type for output parameters list.";

      auto* pname = python::element(pout, MODULE_NAME);  // borrowed
      if (!pname)
        PY_ERROR << "Output module name is required.";
      rt_params_->addEventExporter(
          EventExporterFactory::get().build(python::get<std::string>(pname), python::get<ParametersList>(pout)));
    }

    void PythonHandler::parseExtraParticles(PyObject* pparts) {
      if (!python::is<ParametersList>(pparts))
        PY_ERROR << "Extra particles definition object should be a parameters list.";

      const auto& parts = python::get<ParametersList>(pparts);
      for (const auto& k : parts.keys(true)) {
        auto props = parts.get<ParametersList>(k);
        if (props.has<int>("pdgid"))
          props.set<pdgid_t>("pdgid", props.get<int>("pdgid"));
        const ParticleProperties part(props);
        if (part.mass <= 0. && part.width <= 0.)  // skip aliases
          continue;
        if (!PDG::get().has(part.pdgid) || PDG::get()(part.pdgid) != part) {
          CG_INFO("PythonHandler:particles") << "Adding a new particle with PDG id=" << part.pdgid << " and name \""
                                             << part.name << "\" to the PDG dictionary.";
          PDG::get().define(part);
        }
      }
    }

    ParametersDescription PythonHandler::description() {
      auto desc = Handler::description();
      desc.setDescription("Python 2/3 cards parser");
      desc.add<int>("debugging", 0).setDescription("debugging level");
      desc.add<int>("verbosity", 0).setDescription("verbosity level");
      return desc;
    }
  }  // namespace card
}  // namespace cepgen
typedef cepgen::card::PythonHandler PythonCardHandler;
REGISTER_CARD_HANDLER(".py", PythonCardHandler);
