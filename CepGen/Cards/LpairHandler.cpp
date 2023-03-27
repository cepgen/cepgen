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

#include <fstream>

#include "CepGen/Cards/LpairHandler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"  // for library loading
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Process/Process.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

namespace cepgen {
  namespace card {
    LpairHandler::LpairHandler(const ParametersList& params)
        : Handler(params),
          proc_params_(new ParametersList),
          kin_params_(new ParametersList),
          gen_params_(new ParametersList),
          int_params_(new ParametersList),
          pdg_input_path_("mass_width_2021.mcd") {}

    void LpairHandler::init() {
      //-------------------------------------------------------------------------------------------
      // Process/integration/hadronisation parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<std::string>("PROC", "Process name to simulate", &proc_name_);
      registerIntegratorParameter<std::string>("ITYP", "Integration algorithm", MODULE_NAME);
      registerParameter<std::string>("HADR", "Hadronisation algorithm", &evt_mod_name_);
      registerParameter<std::string>("EVMD", "Events modification algorithms", &evt_mod_name_);
      registerParameter<std::string>("OUTP", "Output module", &out_mod_name_);
      registerParameter<std::string>("OUTF", "Output file name", &out_file_name_);
      registerParameter<std::string>("ADDN", "Additional libraries to load", &addons_list_);

      //-------------------------------------------------------------------------------------------
      // General parameters
      //-------------------------------------------------------------------------------------------

      registerIntegratorParameter<int>("NTRT", "Smoothen the integrand", "treat");
      registerParameter<int>("TIMR", "Enable the time ticker", &timer_);
      registerParameter<int>("IEND", "Generation type", &iend_);
      registerParameter<int>("DEBG", "Debugging verbosity", &log_level_);
      registerParameter<int>("LOGE", "Extended logging", &ext_log_);
      registerIntegratorParameter<int>("NCVG", "Number of function calls", "numFunctionCalls");
      registerIntegratorParameter<int>("ITVG", "Number of integration iterations", "iterations");
      registerIntegratorParameter<int>("SEED", "Random generator seed", "seed");
      registerKinematicsParameter<int>("MODE", "Subprocess' mode", "mode");
      registerGenerationParameter<int>("NTHR", "Number of threads to use for events generation", "numThreads");
      registerGenerationParameter<int>("NCSG", "Number of points to probe", "numPoints");
      registerGenerationParameter<int>("NGEN", "Number of events to generate", "maxgen");
      registerGenerationParameter<int>("NPRN", "Number of events before printout", "printEvery");

      //-------------------------------------------------------------------------------------------
      // Process-specific parameters
      //-------------------------------------------------------------------------------------------

      registerProcessParameter<int>("METH", "Computation method (kT-factorisation)", "method");
      registerProcessParameter<int>("IPOL", "Polarisation states to consider", "polarisationStates");

      //-------------------------------------------------------------------------------------------
      // Process kinematics parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<std::string>("KMRG", "KMR grid interpolation path", &kmr_grid_path_);
      registerParameter<std::string>("MGRD", "MSTW grid interpolation path", &mstw_grid_path_);
      registerParameter<std::string>("PDGI", "Input file for PDG information", &pdg_input_path_);
      registerParameter<int>("PMOD", "Outgoing primary particles' mode", &str_fun_);
      registerParameter<int>("EMOD", "Outgoing primary particles' mode", &str_fun_);
      registerParameter<int>("RTYP", "R-ratio computation type", &sr_type_);
      registerProcessParameter<int>("PAIR", "Outgoing particles' PDG id", "pair");
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

      registerIntegratorParameter<int>("NTREAT", "Smoothen the integrand", "treat");
      registerIntegratorParameter<int>("ITMX", "Number of integration iterations", "iterations");
      registerIntegratorParameter<int>("NCVG", "Number of function calls to perform", "numFunctionCalls");
      registerProcessParameter<int>("METHOD", "Computation method (kT-factorisation)", "method");
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

    Parameters* LpairHandler::parseFile(const std::string& filename, Parameters* params) {
      if (!utils::fileExists(filename))
        throw CG_FATAL("LpairHandler") << "Unable to locate steering card \"" << filename << "\".";
      rt_params_ = params;
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
          if (utils::ltrim(key)[0] == '#')
            continue;
          setParameter(key, value);
          if (describe(key) != kInvalidStr)
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
        pdg::MCDFileParser::parse(pdg_input_path_);
      if (!kmr_grid_path_.empty())
        kmr::GluonGrid::get(ParametersList().set<std::string>("path", kmr_grid_path_));

      //--- build the ticker if required
      if (timer_)
        rt_params_->setTimeKeeper(new utils::TimeKeeper);
      utils::Logger::get().setLevel((utils::Logger::Level)log_level_);
      utils::Logger::get().setExtended(ext_log_);

      //--- parse the structure functions code
      auto sf_params = StructureFunctionsFactory::get().describeParameters(str_fun_).parameters();
      sf_params.set<ParametersList>("sigmaRatio", SigmaRatiosFactory::get().describeParameters(sr_type_).parameters());
      if (str_fun_ == 205 /* MSTWgrid */ && !mstw_grid_path_.empty())
        sf_params.set<std::string>("gridPath", mstw_grid_path_);
      kin_params_->set("structureFunctions", sf_params);
      proc_params_->set("kinematics", *kin_params_);

      //--- parse the process name
      if (!proc_name_.empty() || !proc_params_->empty()) {
        if (!rt_params_->hasProcess() && proc_name_.empty())
          throw CG_FATAL("LpairHandler") << "Process name not specified!";
        if (rt_params_->hasProcess() && rt_params_->process().name() == proc_name_)
          *proc_params_ = ParametersList(rt_params_->process().parameters()) + *proc_params_;
        if (proc_name_ == "pptoff" && lepton_id_ != 0)
          proc_params_->operator[]<int>("pair") = 11 + (lepton_id_ - 1) * 2;
        rt_params_->setProcess(ProcessFactory::get().build(proc_name_, *proc_params_));
      }

      rt_params_->generation().setParameters(*gen_params_);

      rt_params_->par_integrator += *int_params_;

      //--- parse the hadronisation algorithm name
      if (!evt_mod_name_.empty())
        for (const auto& mod : utils::split(evt_mod_name_, ','))
          rt_params_->addModifier(EventModifierFactory::get().build(mod, ParametersList()));

      //--- parse the output module name
      if (!out_mod_name_.empty()) {
        const auto& out_files = utils::split(out_file_name_, ',');
        size_t i = 0;
        for (const auto& mod : utils::split(out_mod_name_, ',')) {
          ParametersList outm;
          if (out_files.size() > i && !out_files.at(i).empty())
            outm.set<std::string>("filename", out_files.at(i));
          rt_params_->addEventExporter(EventExporterFactory::get().build(mod, outm));
          ++i;
        }
      }

      return rt_params_;
    }

    void LpairHandler::write(const std::string& file) const {
      std::map<std::string, std::string> out_map;
      for (const auto& it : p_strings_)
        if (it.second.value && !it.second.value->empty())
          out_map[it.first] = utils::format(
              "%-8s %-20s ! %s\n", it.first.data(), it.second.value->data(), it.second.description.data());
      for (const auto& it : p_ints_)
        if (it.second.value && *it.second.value != kInvalidInt)
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
      rt_params_ = const_cast<Parameters*>(params);
      str_fun_ = rt_params_->kinematics().incomingBeams().structureFunctions().name<int>();
      sr_type_ = rt_params_->kinematics().incomingBeams().structureFunctions().get<int>("sigmaRatio");
      //kmr_grid_path_ = kmr::GluonGrid::get().path();
      //mstw_grid_path_ =
      //pdg_input_path_ =
      iend_ = (int)rt_params_->generation().enabled();
      log_level_ = (int)utils::Logger::get().level();
      ext_log_ = utils::Logger::get().extended();
      proc_name_ = rt_params_->processName();
      *proc_params_ += rt_params_->process().parameters();
      if (proc_params_->has<ParticleProperties>("pair"))
        proc_params_->set<int>("pair", proc_params_->get<ParticleProperties>("pair").pdgid);
      if (proc_name_ == "pptoff" || proc_name_ == "pptoll" /* legacy */)
        lepton_id_ = (rt_params_->process().parameters().get<int>("pair") - 11) / 2. + 1;
      {
        std::vector<std::string> evt_mod;
        std::transform(rt_params_->eventModifiersSequence().begin(),
                       rt_params_->eventModifiersSequence().end(),
                       std::back_inserter(evt_mod),
                       [](const auto& mod) { return mod->name(); });
        evt_mod_name_ = utils::merge(evt_mod, ",");
      }
      {
        std::vector<std::string> out_mod, out_mod_file;
        for (const auto& out : rt_params_->eventExportersSequence()) {
          out_mod.emplace_back(out->name());
          out_mod_file.emplace_back(out->parameters().get<std::string>("filename"));
        }
        out_mod_name_ = utils::merge(out_mod, ",");
        out_file_name_ = utils::merge(out_mod_file, ",");
      }
      timer_ = (rt_params_->timeKeeper() != nullptr);

      *kin_params_ += rt_params_->kinematics().parameters(true);
      *gen_params_ += rt_params_->generation().parameters();
      *int_params_ += rt_params_->par_integrator;
      init();
    }

    void LpairHandler::setParameter(const std::string& key, const std::string& value) {
      // particular case for the double as we cannot rely on casting exceptions
      if (value.find('.') != std::string::npos)
        try {
          set<double>(key, std::stod(value));
          return;
        } catch (const std::logic_error&) {
          for (const auto& let : value)
            if (isalpha(let) && let != 'E' && let != 'e') {
              set<std::string>(key, value);
              return;
            }
          throw CG_FATAL("LpairHandler:setParameter")
              << "Failed to parse a floating-point parameter \"" << key << "\" → \"" << value << "\"!";
        }
      try {
        set<int>(key, std::stoi(value));
      } catch (const std::logic_error&) {
        try {
          set<std::string>(key, value);
        } catch (const std::logic_error&) {
          throw CG_FATAL("LpairHandler:setParameter")
              << "Failed to add the parameter \"" << key << "\" → \"" << value << "\"!";
        }
      }
    }

    std::string LpairHandler::parameter(std::string key) const {
      {
        auto var = get<double>(key);
        if (var != -kInvalidDbl)
          return std::to_string(var);
      }
      {
        auto var = get<int>(key);
        if (var != -kInvalidInt)
          return std::to_string(var);
      }
      return get<std::string>(key);
    }

    std::string LpairHandler::describe(std::string key) const {
      if (p_strings_.count(key))
        return p_strings_.find(key)->second.description;
      if (p_ints_.count(key))
        return p_ints_.find(key)->second.description;
      if (p_doubles_.count(key))
        return p_doubles_.find(key)->second.description;
      return kInvalidStr;
    }

    ParametersDescription LpairHandler::description() {
      auto desc = Handler::description();
      desc.setDescription("LPAIR-like cards parser");
      return desc;
    }
  }  // namespace card
}  // namespace cepgen
typedef cepgen::card::LpairHandler LpairCardHandler;
REGISTER_CARD_HANDLER(".card", LpairCardHandler);
