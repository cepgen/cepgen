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
#include "CepGen/Processes/Process.h"
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

      static ParametersDescription description();

      Parameters* parse(const std::string&, Parameters*) override;

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
    };

    PythonHandler::PythonHandler(const ParametersList& params) : Handler(params) {}

    Parameters* PythonHandler::parse(const std::string& file, Parameters* params) {
      if (!utils::fileExists(file))
        throw CG_FATAL("PythonHandler") << "Unable to locate steering card \"" << file << "\".";
      utils::env::set("PYTHONDONTWRITEBYTECODE", "1");

      rt_params_ = params;
      std::string filename = python::pythonPath(file);
      const size_t fn_len = filename.length() + 1;

      //Py_DebugFlag = 1;
      //Py_VerboseFlag = 1;

      {  // scope of the filename definition
#ifdef PYTHON2
        char* sfilename = new char[fn_len];
        snprintf(sfilename, fn_len, "%s", filename.c_str());
#else
        wchar_t* sfilename = new wchar_t[fn_len];
        swprintf(sfilename, fn_len, L"%s", filename.c_str());
#endif
        if (!sfilename)
          throw CG_FATAL("PythonHandler") << "Invalid filename provided to the Python cards parser!";
        Py_SetProgramName(sfilename);
        delete[] sfilename;
      }

      for (const auto& path : std::vector<std::string>{utils::env::get("CEPGEN_PATH", "."),
                                                       fs::current_path(),
                                                       fs::current_path() / "Cards",
                                                       fs::current_path().parent_path() / "Cards",
                                                       fs::current_path().parent_path().parent_path() / "Cards",
                                                       "/usr/share/CepGen/Cards"})
        utils::env::append("PYTHONPATH", path);

      Py_InitializeEx(1);
      if (!Py_IsInitialized())
        throw CG_FATAL("PythonHandler") << "Failed to initialise the Python cards parser!";

      CG_DEBUG("PythonHandler").log([](auto& log) {
        auto* py_home = Py_GetPythonHome();
#ifdef PYTHON2
        std::string python_path{Py_GetPath()}, python_home{py_home ? py_home : "(not set)"};
#else
        std::wstring python_path{Py_GetPath()}, python_home{py_home ? py_home : L"(not set)"};
#endif
        log << "Initialised the Python cards parser\n\t"
            << "Python version: " << utils::replace_all(std::string{Py_GetVersion()}, "\n", " ") << "\n\t"
            << "Platform: " << Py_GetPlatform() << "\n\t"
            << "Home directory: " << python_home << "\n\t"
            << "Parsed path: " << python_path << ".";
      });

      auto* cfg = PyImport_ImportModule(filename.c_str());  // new
      if (!cfg)
        PY_ERROR << "Failed to import the configuration card '" << filename << "'\n (parsed from '" << file << "').";

      auto parseAttr = [this](PyObject* cfg, const std::string& name, std::function<void(PyObject*)> callback) -> void {
        if (PyObject_HasAttrString(cfg, name.c_str()) != 1)
          return;
        auto pobj = python::ObjectPtr(PyObject_GetAttrString(cfg, name.c_str()));  // new
        if (pobj)
          callback(pobj.get());
      };

      //--- additional libraries to load
      parseAttr(cfg, ADDONS_NAME, [this](PyObject* padd) {
        for (const auto& lib : python::getVector<std::string>(padd))
          loadLibrary(lib);
      });

      //--- timekeeper definition
      // (currently, does not parse the object, just check its presence)
      parseAttr(cfg, TIMER_NAME, [this](PyObject*) { rt_params_->setTimeKeeper(new utils::TimeKeeper); });

      //--- general particles definition
      parseAttr(
          cfg, MCD_NAME, [this](PyObject* ppdg) { pdg::MCDFileParser::parse(python::get<std::string>(ppdg).c_str()); });

      //--- additional particles definition
      parseAttr(cfg, PDGLIST_NAME, [this](PyObject* pextp) { parseExtraParticles(pextp); });

      //--- process definition
      parseAttr(cfg, PROCESS_NAME, [this, &file](PyObject* process) {
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

        rt_params_->par_kinematics += pkin;
        if (proc_params.has<int>("mode"))
          rt_params_->par_kinematics.set<int>("mode", proc_params.get<int>("mode"));

        //--- taming functions
        auto* ptam = python::element(process, "tamingFunctions");  // borrowed
        if (ptam)
          for (const auto& p : python::getVector<ParametersList>(ptam))
            rt_params_->addTamingFunction(utils::FunctionalFactory::get().build("ROOT", p));
      });

      parseAttr(cfg, LOGGER_NAME, [this](PyObject* plog) { parseLogging(plog); });

      //--- hadroniser parameters (legacy)
      parseAttr(cfg, HADR_NAME, [this](PyObject* phad) { parseHadroniser(phad); });
      parseAttr(cfg, EVT_MOD_SEQ_NAME, [this](PyObject* pmod_seq) { parseEventModifiers(pmod_seq); });

      //--- generation parameters
      parseAttr(cfg, INTEGRATOR_NAME, [this](PyObject* pint) { parseIntegrator(pint); });
      parseAttr(cfg, GENERATOR_NAME, [this](PyObject* pgen) { parseGenerator(pgen); });
      parseAttr(cfg, OUTPUT_NAME, [this](PyObject* pout) { parseOutputModules(pout); });

      //--- finalisation
      Py_CLEAR(cfg);

      if (Py_IsInitialized())
        Py_Finalize();

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
