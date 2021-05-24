#include "CepGen/Cards/LpairHandler.h"
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
#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"
#include "CepGen/Utils/Filesystem.h"

#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/GluonGrid.h"

#include <fstream>

namespace cepgen {
  namespace card {
    const int LpairHandler::kInvalid = 99999;

    //----- specialization for LPAIR input cards

    LpairHandler::LpairHandler(const ParametersList& params)
        : Handler(params),
          proc_params_(new ParametersList),
          kin_params_(new ParametersList),
          gen_params_(new ParametersList),
          timer_(false),
          str_fun_(11),
          sr_type_(1),
          lepton_id_(0),
          pdg_input_path_("mass_width_2020.mcd"),
          iend_(1) {}

    Parameters* LpairHandler::parse(const std::string& filename, Parameters* params) {
      if (!utils::fileExists(filename))
        throw CG_FATAL("LpairHandler") << "Unable to locate steering card \"" << filename << "\".";
      params_ = params;
      std::ostringstream os;
      {  //--- file parsing part
        std::ifstream file(filename, std::fstream::in);
        if (!file.is_open())
          throw CG_FATAL("LpairHandler") << "Failed to parse file \"" << filename << "\".";

        init();

        //--- parse all fields
        std::string line, key, value;
        while (getline(file, line)) {
          std::istringstream iss(line);
          iss >> key >> value;
          if (key[0] == '#')  // FIXME need to ensure there is no extra space before!
            continue;
          setParameter(key, value);
          if (describe(key) != "null")
            os << utils::format("\n>> %-8s %-25s (%s)", key.c_str(), parameter(key).c_str(), describe(key).c_str());
        }
        file.close();
      }

      CG_INFO("LpairHandler") << "File '" << filename << "' successfully retrieved!\n\t"
                              << "The following parameters are set:" << os.str() << "\n\t"
                              << "Now parsing the configuration.";

      if (!addons_list_.empty())
        for (const auto& lib : utils::split(addons_list_, ','))
          loadLibrary(lib);

      //--- parse the PDG library
      if (!pdg_input_path_.empty())
        pdg::MCDFileParser::parse(pdg_input_path_.c_str());
      if (!kmr_grid_path_.empty())
        kmr::GluonGrid::get(kmr_grid_path_);

      //--- build the ticker if required
      if (timer_)
        params_->setTimeKeeper(new utils::TimeKeeper);

      //--- parse the process name
      if (!proc_name_.empty() || !proc_params_->empty()) {
        if (!params_->hasProcess() && proc_name_.empty())
          throw CG_FATAL("LpairHandler") << "Process name not specified!";
        if (params_->hasProcess() && params_->process().name() == proc_name_)
          *proc_params_ = ParametersList(params_->process().parameters()) + *proc_params_;
        if (proc_name_ == "pptoff" && lepton_id_ != 0)
          proc_params_->operator[]<int>("pair") = 11 + (lepton_id_ - 1) * 2;
        params_->setProcess(proc::ProcessesFactory::get().build(proc_name_, *proc_params_));
      }

      params_->kinematics = Kinematics(*kin_params_);
      params_->generation() = Parameters::Generation(*gen_params_);

      //--- parse the structure functions code
      if (str_fun_ == (int)strfun::Type::MSTWgrid && !mstw_grid_path_.empty())
        params_->kinematics.incoming_beams.setStructureFunctions(strfun::StructureFunctionsFactory::get().build(
            str_fun_, ParametersList().set<std::string>("gridPath", mstw_grid_path_)));
      else
        params_->kinematics.incoming_beams.setStructureFunctions(str_fun_, sr_type_);

      //--- parse the hadronisation algorithm name
      if (!evt_mod_name_.empty())
        for (const auto& mod : utils::split(evt_mod_name_, ','))
          params_->addModifier(EventModifierFactory::get().build(mod, ParametersList()));

      //--- parse the output module name
      if (!out_mod_name_.empty()) {
        const auto& out_files = utils::split(out_file_name_, ',');
        size_t i = 0;
        for (const auto& mod : utils::split(out_mod_name_, ',')) {
          ParametersList outm;
          if (out_files.size() > i && !out_files.at(i).empty())
            outm.set<std::string>("filename", out_files.at(i));
          params_->addOutputModule(io::ExportModuleFactory::get().build(mod, outm));
          ++i;
        }
      }

      return params_;
    }

    void LpairHandler::init() {
      //-------------------------------------------------------------------------------------------
      // Process/integration/hadronisation parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<std::string>("PROC", "Process name to simulate", &proc_name_);
      registerParameter<std::string>(
          "ITYP", "Integration algorithm", &params_->integrator->operator[]<std::string>(ParametersList::MODULE_NAME));
      registerParameter<std::string>("HADR", "Hadronisation algorithm", &evt_mod_name_);
      registerParameter<std::string>("EVMD", "Events modification algorithms", &evt_mod_name_);
      registerParameter<std::string>("OUTP", "Output module", &out_mod_name_);
      registerParameter<std::string>("OUTF", "Output file name", &out_file_name_);
      registerParameter<std::string>("ADDN", "Additional libraries to load", &addons_list_);

      //-------------------------------------------------------------------------------------------
      // General parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<int>("NTRT", "Smoothen the integrand", (int*)&params_->integrator->operator[]<bool>("treat"));
      registerParameter<int>("TIMR", "Enable the time ticker", &timer_);
      registerParameter<int>("IEND", "Generation type", &iend_);
      registerParameter<int>("DEBG", "Debugging verbosity", (int*)&utils::Logger::get().level);
      registerParameter<int>(
          "NCVG", "Number of function calls", (int*)&params_->integrator->operator[]<int>("numFunctionCalls"));
      registerParameter<int>(
          "ITVG", "Number of integration iterations", (int*)&params_->integrator->operator[]<int>("iterations"));
      registerParameter<int>("SEED", "Random generator seed", (int*)&params_->integrator->operator[]<int>("seed"));
      registerParameter<int>("MODE", "Subprocess' mode", &kin_params_->operator[]<int>("mode"));
      registerParameter<int>(
          "NTHR", "Number of threads to use for events generation", &gen_params_->operator[]<int>("numThreads"));
      registerParameter<int>("NCSG", "Number of points to probe", &gen_params_->operator[]<int>("numPoints"));
      registerParameter<int>("NGEN", "Number of events to generate", &gen_params_->operator[]<int>("maxgen"));
      registerParameter<int>("NPRN", "Number of events before printout", &gen_params_->operator[]<int>("printEvery"));

      //-------------------------------------------------------------------------------------------
      // Process-specific parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<int>("METH", "Computation method (kT-factorisation)", &proc_params_->operator[]<int>("method"));
      registerParameter<int>(
          "IPOL", "Polarisation states to consider", &proc_params_->operator[]<int>("polarisationStates"));

      //-------------------------------------------------------------------------------------------
      // Process kinematics parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<std::string>("KMRG", "KMR grid interpolation path", &kmr_grid_path_);
      registerParameter<std::string>("MGRD", "MSTW grid interpolation path", &mstw_grid_path_);
      registerParameter<std::string>("PDGI", "Input file for PDG information", &pdg_input_path_);
      registerParameter<int>("PMOD", "Outgoing primary particles' mode", &str_fun_);
      registerParameter<int>("EMOD", "Outgoing primary particles' mode", &str_fun_);
      registerParameter<int>("RTYP", "R-ratio computation type", &sr_type_);
      registerParameter<int>("PAIR", "Outgoing particles' PDG id", (int*)&proc_params_->operator[]<int>("pair"));
      registerKinematicsParameter<std::string>("FFAC", "Form factors for the incoming beams", "formFactors");
      registerKinematicsParameter<int>("INA1", "Heavy ion atomic weight (1st incoming beam)", "beam1A");
      registerKinematicsParameter<int>("INZ1", "Heavy ion atomic number (1st incoming beam)", "beam1Z");
      registerKinematicsParameter<int>("INA2", "Heavy ion atomic weight (2nd incoming beam)", "beam2A");
      registerKinematicsParameter<int>("INZ2", "Heavy ion atomic number (2nd incoming beam)", "beam2Z");
      registerKinematicsParameter<double>("INP1", "Momentum (1st primary particle)", "beam1pz");
      registerKinematicsParameter<double>("INP2", "Momentum (2nd primary particle)", "beam2pz");
      registerKinematicsParameter<double>("INPP", "Momentum (1st primary particle)", "beam1pz");
      registerKinematicsParameter<double>("INPE", "Momentum (2nd primary particle)", "beam2pz");
      registerKinematicsParameter<double>(
          "PTCT", "Minimal transverse momentum (single central outgoing particle)", "ptmin");
      registerKinematicsParameter<double>(
          "PTMX", "Maximal transverse momentum (single central outgoing particle)", "ptmax");
      registerKinematicsParameter<double>("MSCT", "Minimal central system mass", "invmassmin");
      registerKinematicsParameter<double>("MSMX", "Maximal central system mass", "invmassmax");
      registerKinematicsParameter<double>("ECUT", "Minimal energy (single central outgoing particle)", "energysummin");
      registerKinematicsParameter<double>("ETMN", "Minimal pseudo-rapidity (central outgoing particles)", "etamin");
      registerKinematicsParameter<double>("ETMX", "Maximal pseudo-rapidity (central outgoing particles)", "etamax");
      registerKinematicsParameter<double>("YMIN", "Minimal rapidity (central outgoing particles)", "rapiditymin");
      registerKinematicsParameter<double>("YMAX", "Maximal rapidity (central outgoing particles)", "rapiditymax");
      registerKinematicsParameter<double>(
          "PDMN", "Minimal transverse momentum difference (central outgoing particles)", "ptdiffmin");
      registerKinematicsParameter<double>(
          "PDMX", "Maximal transverse momentum difference (central outgoing particles)", "ptdiffmax");
      registerKinematicsParameter<double>("Q2MN", "Minimal Q^2 = -q^2 (exchanged parton)", "q2min");
      registerKinematicsParameter<double>("Q2MX", "Maximal Q^2 = -q^2 (exchanged parton)", "q2max");
      registerKinematicsParameter<double>("QTMN", "Minimal Q_T (exchanged parton)", "qtmin");
      registerKinematicsParameter<double>("QTMX", "Maximal Q_T (exchanged parton)", "qtmax");
      registerKinematicsParameter<double>("MXMN", "Minimal invariant mass of proton remnants", "mxmin");
      registerKinematicsParameter<double>("MXMX", "Maximal invariant mass of proton remnants", "mxmax");
      registerKinematicsParameter<double>("XIMN", "Minimal fractional momentum loss of outgoing proton (xi)", "ximin");
      registerKinematicsParameter<double>("XIMX", "Maximal fractional momentum loss of outgoing proton (xi)", "ximax");
      registerKinematicsParameter<double>("YJMN", "Minimal remnant jet rapidity", "yjmin");
      registerKinematicsParameter<double>("YJMX", "Maximal remnant jet rapidity", "yjmax");

      //-------------------------------------------------------------------------------------------
      // PPtoLL cards backward compatibility
      //-------------------------------------------------------------------------------------------

      registerParameter<int>("NTREAT", "Smoothen the integrand", (int*)&params_->integrator->operator[]<bool>("treat"));
      registerParameter<int>(
          "ITMX", "Number of integration iterations", (int*)&params_->integrator->operator[]<int>("iterations"));
      registerParameter<int>("NCVG", "Number of points to probe", (int*)&params_->generation().num_points);
      registerParameter<int>(
          "METHOD", "Computation method (kT-factorisation)", &proc_params_->operator[]<int>("method"));
      registerParameter<int>("LEPTON", "Outgoing leptons' flavour", &lepton_id_);
      registerKinematicsParameter<double>(
          "PTMIN", "Minimal transverse momentum (single central outgoing particle)", "ptmin");
      registerKinematicsParameter<double>(
          "PTMAX", "Maximal transverse momentum (single central outgoing particle)", "ptmax");
      registerKinematicsParameter<double>("Q1TMIN", "Minimal Q_T (exchanged parton)", "qtmin");
      registerKinematicsParameter<double>("Q1TMAX", "Maximal Q_T (exchanged parton)", "qtmax");
      registerKinematicsParameter<double>("Q2TMIN", "Minimal Q_T (exchanged parton)", "qtmin");
      registerKinematicsParameter<double>("Q2TMAX", "Maximal Q_T (exchanged parton)", "qtmax");
      registerKinematicsParameter<double>("MXMIN", "Minimal invariant mass of proton remnants", "mxmin");
      registerKinematicsParameter<double>("MXMAX", "Maximal invariant mass of proton remnants", "mxmax");
    }

    void LpairHandler::write(const std::string& file) const {
      std::map<std::string, std::string> out_map;
      for (const auto& it : p_strings_)
        if (it.second.value && !it.second.value->empty())
          out_map[it.first] = utils::format(
              "%-8s %-20s ! %s\n", it.first.data(), it.second.value->data(), it.second.description.data());
      for (const auto& it : p_ints_)
        if (it.second.value && *it.second.value != kInvalid)
          out_map[it.first] =
              utils::format("%-8s %-20d ! %s\n", it.first.data(), *it.second.value, it.second.description.data());
      for (const auto& it : p_doubles_)
        if (it.second.value && *it.second.value != Limits::INVALID)
          out_map[it.first] =
              utils::format("%-8s %-20e ! %s\n", it.first.data(), *it.second.value, it.second.description.data());

      std::ofstream f(file, std::fstream::out | std::fstream::trunc);
      if (!f.is_open())
        throw CG_ERROR("LpairHandler") << "Failed to open file \"" << file << "%s\" for writing.";
      for (const auto& ln : out_map)
        f << ln.second;
      f.close();
    }

    void LpairHandler::pack(const Parameters* params) {
      params_ = const_cast<Parameters*>(params);
      str_fun_ = params_->kinematics.incoming_beams.structureFunctions()->name();
      if (params_->kinematics.incoming_beams.structureFunctions() &&
          params_->kinematics.incoming_beams.structureFunctions()->sigmaRatio())
        sr_type_ = params_->kinematics.incoming_beams.structureFunctions()->sigmaRatio()->name();
      //kmr_grid_path_ =
      //mstw_grid_path_ =
      //pdg_input_path_ =
      iend_ = (int)params_->generation().enabled();
      proc_name_ = params_->processName();
      *proc_params_ += params_->process().parameters();
      if (proc_params_->has<ParticleProperties>("pair"))
        proc_params_->set<int>("pair", proc_params_->get<ParticleProperties>("pair").pdgid);
      if (proc_name_ == "pptoff")
        lepton_id_ = (params_->process().parameters().get<int>("pair") - 11) / 2. + 1;
      {
        std::vector<std::string> evt_mod;
        for (const auto& mod : params_->eventModifiersSequence())
          evt_mod.emplace_back(mod->name());
        evt_mod_name_ = utils::merge(evt_mod, ",");
      }
      {
        std::vector<std::string> out_mod, out_mod_file;
        for (const auto& out : params_->outputModulesSequence()) {
          out_mod.emplace_back(out->name());
          out_mod_file.emplace_back(out->parameters().get<std::string>("filename"));
        }
        out_mod_name_ = utils::merge(out_mod, ",");
        out_file_name_ = utils::merge(out_mod_file, ",");
      }
      timer_ = (params_->timeKeeper() != nullptr);

      *kin_params_ += params_->kinematics.parameters();
      *gen_params_ += params_->generation().parameters();
      init();
    }

    void LpairHandler::setParameter(const std::string& key, const std::string& value) {
      // particular case for the double as we cannot rely on casting exceptions
      if (value.find('.') != std::string::npos)
        try {
          setValue<double>(key.c_str(), std::stod(value));
          return;
        } catch (const std::logic_error&) {
          for (const auto& let : value)
            if (isalpha(let) && let != 'E' && let != 'e') {
              setValue<std::string>(key.c_str(), value);
              return;
            }
          throw CG_FATAL("LpairHandler:setParameter")
              << "Failed to parse a floating-point parameter \"" << key << "\" → \"" << value << "\"!";
        }
      try {
        setValue<int>(key.c_str(), std::stoi(value));
      } catch (const std::logic_error&) {
        try {
          setValue<std::string>(key.c_str(), value);
        } catch (const std::logic_error&) {
          throw CG_FATAL("LpairHandler:setParameter")
              << "Failed to add the parameter \"" << key << "\" → \"" << value << "\"!";
        }
      }
    }

    std::string LpairHandler::parameter(std::string key) const {
      {
        auto var = getValue<double>(key.c_str());
        if (var != -999.)
          return std::to_string(var);
      }
      {
        auto var = getValue<int>(key.c_str());
        if (var != -999999)
          return std::to_string(var);
      }
      return getValue<std::string>(key.c_str());
    }

    std::string LpairHandler::describe(std::string key) const {
      if (p_strings_.count(key))
        return p_strings_.find(key)->second.description;
      if (p_ints_.count(key))
        return p_ints_.find(key)->second.description;
      if (p_doubles_.count(key))
        return p_doubles_.find(key)->second.description;
      return "null";
    }
  }  // namespace card
}  // namespace cepgen

REGISTER_CARD_HANDLER(".card", LpairHandler)
