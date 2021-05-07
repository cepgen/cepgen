#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Generator.h"  // for library loading
#include "CepGen/Parameters.h"

#include "CepGen/Modules/CardsHandlerFactory.h"

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Modules/EventModifierFactory.h"

#include "CepGen/Core/ExportModule.h"
#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Processes/Process.h"
#include "CepGen/Modules/ProcessesFactory.h"

#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include "CepGen/Utils/Functional.h"
#include "CepGen/Modules/FunctionalFactory.h"

#include "CepGen/Integration/Integrator.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/HeavyIon.h"

#include "CepGen/Utils/TimeKeeper.h"
#include "CepGen/Utils/String.h"

#include <algorithm>

#if PY_MAJOR_VERSION < 3
#define PYTHON2
#endif

namespace cepgen {
  namespace card {
    PythonHandler::PythonHandler(const ParametersList& params) : Handler(params) {}

    Parameters* PythonHandler::parse(const std::string& file, Parameters* params) {
      setenv(
          "PYTHONPATH", (utils::environ("CEPGEN_PATH", ".") + ":.:Cards:../Cards:/usr/share/CepGen/Cards").c_str(), 1);
      setenv("PYTHONDONTWRITEBYTECODE", "1", 1);
      CG_DEBUG("PythonHandler") << "Python PATH: \"" << utils::environ("PYTHONPATH") << "\".";

      params_ = params;
      std::string filename = pythonPath(file);
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

      Py_InitializeEx(1);

      if (!Py_IsInitialized())
        throw CG_FATAL("PythonHandler") << "Failed to initialise the Python cards parser!";

      CG_DEBUG("PythonHandler") << "Initialised the Python cards parser\n\t"
                                << "Python version: " << Py_GetVersion() << "\n\t"
                                << "Platform: " << Py_GetPlatform() << ".";

      PyObject* cfg = PyImport_ImportModule(filename.c_str());  // new
      if (!cfg)
        throwPythonError("Failed to import the configuration card '" + file + "'\n (parsed as '" + filename + "')");

      //--- additional libraries to load
      if (PyObject_HasAttrString(cfg, ADDONS_NAME) == 1) {
        PyObject* padd = PyObject_GetAttrString(cfg, ADDONS_NAME);  // new
        if (padd) {
          for (const auto& lib : getVector<std::string>(padd))
            loadLibrary(lib);
          Py_CLEAR(padd);
        }
      }

      //--- timekeeper definition
      if (PyObject_HasAttrString(cfg, TIMER_NAME) == 1) {
        PyObject* ptim = PyObject_GetAttrString(cfg, TIMER_NAME);  // new
        if (ptim) {
          params_->setTimeKeeper(new utils::TimeKeeper);
          Py_CLEAR(ptim);
        }
      }

      //--- general particles definition
      if (PyObject_HasAttrString(cfg, MCD_NAME) == 1) {
        PyObject* ppdg = PyObject_GetAttrString(cfg, MCD_NAME);  // new
        if (ppdg) {
          pdg::MCDFileParser::parse(get<std::string>(ppdg).c_str());
          Py_CLEAR(ppdg);
        }
      }

      //--- additional particles definition
      if (PyObject_HasAttrString(cfg, PDGLIST_NAME) == 1) {
        PyObject* pextp = PyObject_GetAttrString(cfg, PDGLIST_NAME);  // new
        if (pextp) {
          parseExtraParticles(pextp);
          Py_CLEAR(pextp);
        }
      }

      //--- process definition
      PyObject* process = nullptr;
      if (PyObject_HasAttrString(cfg, PROCESS_NAME) != 1 ||
          !(process = PyObject_GetAttrString(cfg, PROCESS_NAME)))  // new
        throwPythonError("Failed to extract a '" + std::string(PROCESS_NAME) +
                         "' keyword from the configuration card '" + file + "'!");

      //--- list of process-specific parameters
      ParametersList proc_params;
      fillParameter(process, "processParameters", proc_params);

      //--- type of process to consider
      PyObject* pproc_name = element(process, ParametersList::MODULE_NAME);  // borrowed
      if (!pproc_name)
        throwPythonError("Failed to extract the process name from the configuration card '" + file + "'!");

      //--- process mode
      params_->setProcess(proc::ProcessesFactory::get().build(get<std::string>(pproc_name), proc_params));

      //--- process kinematics
      ParametersList pkin;
      PyObject* pin_kinematics = element(process, "inKinematics");  // borrowed
      if (pin_kinematics)
        pkin += get<ParametersList>(pin_kinematics);

      PyObject* pout_kinematics = element(process, "outKinematics");  // borrowed
      if (pout_kinematics)
        pkin += get<ParametersList>(pout_kinematics);

      params_->kinematics = Kinematics(pkin);
      if (proc_params.has<int>("mode"))
        params_->kinematics.incoming_beams.setMode((mode::Kinematics)proc_params.get<int>("mode"));

      //--- taming functions
      PyObject* ptam = element(process, "tamingFunctions");  // borrowed
      if (ptam)
        for (const auto& p : getVector<ParametersList>(ptam))
          params_->addTamingFunction(utils::FunctionalFactory::get().build("ROOT", p));

      Py_CLEAR(process);

      if (PyObject_HasAttrString(cfg, LOGGER_NAME) == 1) {
        PyObject* plog = PyObject_GetAttrString(cfg, LOGGER_NAME);  // new
        if (plog) {
          parseLogging(plog);
          Py_CLEAR(plog);
        }
      }

      //--- hadroniser parameters (legacy)
      if (PyObject_HasAttrString(cfg, HADR_NAME) == 1) {
        PyObject* phad = PyObject_GetAttrString(cfg, HADR_NAME);  // new
        if (phad) {
          parseHadroniser(phad);
          Py_CLEAR(phad);
        }
      }

      if (PyObject_HasAttrString(cfg, EVT_MOD_SEQ_NAME) == 1) {
        PyObject* pmod_seq = PyObject_GetAttrString(cfg, EVT_MOD_SEQ_NAME);  // new
        if (pmod_seq) {
          parseEventModifiers(pmod_seq);
          Py_CLEAR(pmod_seq);
        }
      }

      //--- generation parameters
      if (PyObject_HasAttrString(cfg, INTEGRATOR_NAME) == 1) {
        PyObject* pint = PyObject_GetAttrString(cfg, INTEGRATOR_NAME);  // new
        if (pint) {
          parseIntegrator(pint);
          Py_CLEAR(pint);
        }
      }

      if (PyObject_HasAttrString(cfg, GENERATOR_NAME) == 1) {
        PyObject* pgen = PyObject_GetAttrString(cfg, GENERATOR_NAME);  // new
        if (pgen) {
          parseGenerator(pgen);
          Py_CLEAR(pgen);
        }
      }

      if (PyObject_HasAttrString(cfg, OUTPUT_NAME) == 1) {
        PyObject* pout = PyObject_GetAttrString(cfg, OUTPUT_NAME);  // new
        if (pout) {
          if (isVector<ParametersList>(pout))
            parseOutputModules(pout);
          else
            parseOutputModule(pout);
          Py_CLEAR(pout);
        }
      }

      //--- finalisation
      Py_CLEAR(cfg);

      if (Py_IsInitialized())
        Py_Finalize();

      return params_;
    }

    void PythonHandler::parseLogging(PyObject* log) {
      int log_level = 0;
      fillParameter(log, "level", log_level);
      utils::Logger::get().level = (utils::Logger::Level)log_level;
      std::vector<std::string> enabled_modules;
      fillParameter(log, "enabledModules", enabled_modules);
      for (const auto& mod : enabled_modules)
        utils::Logger::get().addExceptionRule(mod);
    }

    void PythonHandler::parseIntegrator(PyObject* integr) {
      if (!PyDict_Check(integr))
        throwPythonError("Integrator object should be a dictionary!");
      *params_->integrator = get<ParametersList>(integr);
    }

    void PythonHandler::parseGenerator(PyObject* gen) {
      if (!PyDict_Check(gen))
        throwPythonError("Generation information object should be a dictionary!");
      params_->generation().enabled = true;
      fillParameter(gen, "numEvents", params_->generation().maxgen);
      fillParameter(gen, "printEvery", params_->generation().gen_print_every);
      fillParameter(gen, "numThreads", params_->generation().num_threads);
      fillParameter(gen, "numPoints", params_->generation().num_points);
    }

    void PythonHandler::parseEventModifiers(PyObject* mod) {
      if (!PyList_Check(mod))
        throwPythonError("Event modification definition object should be a list/Sequence!");

      for (Py_ssize_t i = 0; i < PyList_Size(mod); ++i)
        parseHadroniser(PyList_GetItem(mod, i));
    }

    void PythonHandler::parseHadroniser(PyObject* mod) {
      if (!PyDict_Check(mod))
        throwPythonError("Event modification definition object should be a dictionary!");

      PyObject* pname = element(mod, ParametersList::MODULE_NAME);  // borrowed
      if (!pname)
        throwPythonError("Event modification algorithm name is required!");
      std::string mod_name = get<std::string>(pname);

      params_->addModifier(EventModifierFactory::get().build(mod_name, get<ParametersList>(mod)));

      auto h = params_->eventModifiersSequence().rbegin()->get();
      {  //--- before calling the init() method
        std::vector<std::string> config;
        fillParameter(mod, "preConfiguration", config);
        h->readStrings(config);
      }
      h->init();
      {  //--- after init() has been called
        std::vector<std::string> config;
        fillParameter(mod, "processConfiguration", config);
        for (const auto& block : config) {
          std::vector<std::string> config_blk;
          fillParameter(mod, block.c_str(), config_blk);
          h->readStrings(config_blk);
        }
      }
    }

    void PythonHandler::parseOutputModules(PyObject* mod) {
      if (!PyList_Check(mod))
        throwPythonError("Output modules definition object should be a list/Sequence!");

      for (Py_ssize_t i = 0; i < PyList_Size(mod); ++i)
        parseOutputModule(PyList_GetItem(mod, i));
    }

    void PythonHandler::parseOutputModule(PyObject* pout) {
      if (!is<ParametersList>(pout))
        throwPythonError("Invalid type for output parameters list!");

      PyObject* pname = element(pout, ParametersList::MODULE_NAME);  // borrowed
      if (!pname)
        throwPythonError("Output module name is required!");
      params_->addOutputModule(
          io::ExportModuleFactory::get().build(get<std::string>(pname), get<ParametersList>(pout)));
    }

    void PythonHandler::parseExtraParticles(PyObject* pparts) {
      if (!is<ParametersList>(pparts))
        throwPythonError("Extra particles definition object should be a parameters list!");

      const auto& parts = get<ParametersList>(pparts);
      for (const auto& k : parts.keys(true)) {
        const auto& part = parts.get<ParticleProperties>(k);
        if (part.pdgid == 0 || part.mass < 0.)
          continue;
        CG_DEBUG("PythonHandler:particles")
            << "Adding a new particle with name \"" << part.name << "\" to the PDG dictionary.";
        PDG::get().define(part);
      }
    }
  }  // namespace card
}  // namespace cepgen

REGISTER_CARD_HANDLER("py", PythonHandler)
