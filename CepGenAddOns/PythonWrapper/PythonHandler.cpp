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

// clang-format off
#include "CepGenAddOns/PythonWrapper/PythonError.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"
// clang-format on

#include <algorithm>

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Generator.h"  // for library loading
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

#define Py_DEBUG

namespace cepgen {
  namespace card {
    /// CepGen Python configuration cards reader/writer
    class PythonHandler final : public Handler {
    public:
      /// Read a standard configuration card
      explicit PythonHandler(const ParametersList&);
      Parameters* parse(const std::string&, Parameters*) override;

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

      void parseLogging(PyObject*);
      void parseIntegrator(PyObject*);
      void parseGenerator(PyObject*);
      void parseHadroniser(PyObject*);
      void parseEventModifiers(PyObject*);
      void parseOutputModule(PyObject*);
      void parseOutputModules(PyObject*);
      void parseExtraParticles(PyObject*);
      python::Environment env_;
    };

    PythonHandler::PythonHandler(const ParametersList& params) : Handler(params) {
      //Py_DebugFlag = 1;
      //Py_VerboseFlag = 1;
    }

    Parameters* PythonHandler::parse(const std::string& file, Parameters* params) {
      if (!utils::fileExists(file))
        throw CG_FATAL("PythonHandler") << "Unable to locate steering card \"" << file << "\".";

      rt_params_ = params;
      std::string filename = python::pythonPath(file);
      python::setProgramName(filename);

      CG_DEBUG("PythonHandler").log([](auto& log) {
        log << "Initialised the Python cards parser.";
        for (const auto& ln : python::info())
          log << "\n\t" << ln;
      });

      auto cfg = python::importModule(filename);  // new
      if (!cfg)
        throw PY_ERROR << "Failed to import the configuration card '" << filename << "'\n (parsed from '" << file
                       << "').";

      auto parseAttr = [this, &cfg](const std::string& name, std::function<void(PyObject*)> callback) -> void {
        auto pobj = python::getAttribute(cfg, name);
        if (pobj)
          callback(pobj.get());
      };

      //--- additional libraries to load
      parseAttr(ADDONS_NAME, [this](PyObject* padd) {
        for (const auto& lib : python::getVector<std::string>(padd))
          loadLibrary(lib);
      });

      //--- timekeeper definition
      // (currently, does not parse the object, just check its presence)
      parseAttr(TIMER_NAME, [this](PyObject*) { rt_params_->setTimeKeeper(new utils::TimeKeeper); });

      //--- general particles definition
      parseAttr(MCD_NAME,
                [this](PyObject* ppdg) { pdg::MCDFileParser::parse(python::get<std::string>(ppdg).c_str()); });

      //--- additional particles definition
      parseAttr(PDGLIST_NAME, [this](PyObject* pextp) { parseExtraParticles(pextp); });

      //--- process definition
      parseAttr(PROCESS_NAME, [this, &file](PyObject* process) {
        //--- list of process-specific parameters
        ParametersList proc_params;
        python::fillParameter(process, "processParameters", proc_params);

        //--- type of process to consider
        auto* pproc_name = python::element(process, ParametersList::MODULE_NAME);  // borrowed
        if (!pproc_name)
          PY_ERROR << "Failed to extract the process name from the configuration card '" << file << "'.";

        //--- process mode
        rt_params_->setProcess(proc::ProcessFactory::get().build(python::get<std::string>(pproc_name), proc_params));

        //--- process kinematics
        ParametersList pkin;
        auto* pin_kinematics = python::element(process, "inKinematics");  // borrowed
        if (pin_kinematics)
          pkin += python::get<ParametersList>(pin_kinematics);

        auto* pout_kinematics = python::element(process, "outKinematics");  // borrowed
        if (pout_kinematics)
          pkin += python::get<ParametersList>(pout_kinematics);

        if (proc_params.has<int>("mode"))
          pkin.set<int>("mode", proc_params.get<int>("mode"));
        rt_params_->process().setKinematics(Kinematics(pkin));

        //--- taming functions
        auto* ptam = python::element(process, "tamingFunctions");  // borrowed
        if (ptam)
          for (const auto& p : python::getVector<ParametersList>(ptam))
            rt_params_->addTamingFunction(utils::FunctionalFactory::get().build("ROOT", p));
      });

      parseAttr(LOGGER_NAME, [this](PyObject* plog) { parseLogging(plog); });

      //--- hadroniser parameters (legacy)
      parseAttr(HADR_NAME, [this](PyObject* phad) { parseHadroniser(phad); });
      parseAttr(EVT_MOD_SEQ_NAME, [this](PyObject* pmod_seq) { parseEventModifiers(pmod_seq); });

      //--- generation parameters
      parseAttr(INTEGRATOR_NAME, [this](PyObject* pint) { parseIntegrator(pint); });
      parseAttr(GENERATOR_NAME, [this](PyObject* pgen) { parseGenerator(pgen); });
      parseAttr(OUTPUT_NAME, [this](PyObject* pout) { parseOutputModules(pout); });

      return rt_params_;
    }

    void PythonHandler::parseLogging(PyObject* log) {
      int log_level{(int)utils::Logger::get().level};
      python::fillParameter(log, "level", log_level);
      utils::Logger::get().level = (utils::Logger::Level)log_level;
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

      auto* pname = python::element(mod, ParametersList::MODULE_NAME);  // borrowed
      if (!pname)
        PY_ERROR << "Event modification algorithm name is required.";
      std::string mod_name = python::get<std::string>(pname);

      rt_params_->addModifier(EventModifierFactory::get().build(mod_name, python::get<ParametersList>(mod)));

      auto h = rt_params_->eventModifiersSequence().rbegin()->get();
      {  //--- before calling the init() method
        std::vector<std::string> config;
        python::fillParameter(mod, "preConfiguration", config);
        h->readStrings(config);
      }
      h->init();
      {  //--- after init() has been called
        std::vector<std::string> config;
        python::fillParameter(mod, "processConfiguration", config);
        for (const auto& block : config) {
          std::vector<std::string> config_blk;
          python::fillParameter(mod, block.c_str(), config_blk);
          h->readStrings(config_blk);
        }
      }
    }

    void PythonHandler::parseOutputModules(PyObject* mod) {
      if (!python::isVector<ParametersList>(mod))
        PY_ERROR << "Output modules definition object should be a list/Sequence.";

      for (Py_ssize_t i = 0; i < PyList_Size(mod); ++i)
        parseOutputModule(PyList_GetItem(mod, i));
    }

    void PythonHandler::parseOutputModule(PyObject* pout) {
      if (!python::is<ParametersList>(pout))
        PY_ERROR << "Invalid type for output parameters list.";

      auto* pname = python::element(pout, ParametersList::MODULE_NAME);  // borrowed
      if (!pname)
        PY_ERROR << "Output module name is required.";
      rt_params_->addOutputModule(
          io::ExportModuleFactory::get().build(python::get<std::string>(pname), python::get<ParametersList>(pout)));
    }

    void PythonHandler::parseExtraParticles(PyObject* pparts) {
      if (!python::is<ParametersList>(pparts))
        PY_ERROR << "Extra particles definition object should be a parameters list.";

      const auto& parts = python::get<ParametersList>(pparts);
      for (const auto& k : parts.keys(true)) {
        const auto& part = parts.get<ParticleProperties>(k);
        if (part.pdgid == 0 || part.mass < 0.)
          continue;
        CG_DEBUG("PythonHandler:particles")
            << "Adding a new particle with name \"" << part.name << "\" to the PDG dictionary.";
        PDG::get().define(part);
      }
    }

    ParametersDescription PythonHandler::description() {
      auto desc = Handler::description();
      desc.setDescription("Python 2/3 cards parser");
      return desc;
    }
  }  // namespace card
}  // namespace cepgen

REGISTER_CARD_HANDLER(".py", PythonHandler)
