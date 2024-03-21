/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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
#include <memory>

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
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

using std::string;

namespace cepgen {
  class ParametersList;
  namespace card {
    /// LPAIR-like steering cards parser and writer
    class LpairHandler final : public Handler {
    public:
      /// Read a LPAIR steering card
      explicit LpairHandler(const ParametersList& params)
          : Handler(params),
            proc_params_(new ParametersList),
            kin_params_(new ParametersList),
            gen_params_(new ParametersList),
            int_params_(new ParametersList),
            pdg_input_path_("mass_width_2023.txt") {}

      static ParametersDescription description() {
        auto desc = Handler::description();
        desc.setDescription("LPAIR-like cards parser");
        return desc;
      }

      LpairHandler& parseFile(const std::string& filename) override {
        if (!utils::fileExists(filename))
          throw CG_FATAL("LpairHandler") << "Unable to locate steering card \"" << filename << "\".";
        std::ostringstream os;
        init();
        for (auto line : utils::split(utils::readFile(filename), '\n')) {  // file parsing part
          if (line = utils::trim(line); line.empty() || line[0] == '#')    // skip comments
            continue;
          const auto fields = utils::split(line, ' ', true);  // parse all fields
          if (fields.size() < 2)                              // break on invalid lines
            throw CG_FATAL("LpairHandler") << "Invalid line read from steering card: '" << line << "'.";
          const auto key = fields.at(0), value = fields.at(1);
          setParameter(key, value);
          if (const auto descr = describe(key); descr != kInvalidStr)
            os << utils::format("\n>> %-8s %-25s (%s)", key.data(), parameter(key).data(), descr.data());
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
          runParameters()->setTimeKeeper(new utils::TimeKeeper);
        utils::Logger::get().setLevel((utils::Logger::Level)log_level_);
        utils::Logger::get().setExtended(ext_log_);

        //--- parse the structure functions code
        auto sf_params = StructureFunctionsFactory::get().describeParameters(str_fun_).parameters();
        sf_params.set<ParametersList>("sigmaRatio",
                                      SigmaRatiosFactory::get().describeParameters(sr_type_).parameters());
        if (str_fun_ == 205 /* MSTWgrid */ && !mstw_grid_path_.empty())
          sf_params.set<std::string>("gridPath", mstw_grid_path_);
        kin_params_->set("structureFunctions", sf_params);
        proc_params_->set("kinematics", *kin_params_);

        //--- parse the process name
        if (!proc_name_.empty() || !proc_params_->empty()) {
          if (!runParameters()->hasProcess() && proc_name_.empty())
            throw CG_FATAL("LpairHandler") << "Process name not specified!";
          if (runParameters()->hasProcess() && runParameters()->process().name() == proc_name_)
            *proc_params_ = ParametersList(runParameters()->process().parameters()) + *proc_params_;
          if (proc_name_ == "pptoff" && lepton_id_ != 0)
            proc_params_->operator[]<int>("pair") = 11 + (lepton_id_ - 1) * 2;
          runParameters()->setProcess(ProcessFactory::get().build(proc_name_, *proc_params_));
        }

        runParameters()->integrator() += *int_params_;
        runParameters()->generation().setParameters(*gen_params_);

        //--- parse the hadronisation algorithm name
        if (!evt_mod_name_.empty())
          for (const auto& mod : utils::split(evt_mod_name_, ','))
            runParameters()->addModifier(EventModifierFactory::get().build(mod, ParametersList()));

        //--- parse the output module name
        if (!out_mod_name_.empty()) {
          const auto& out_files = utils::split(out_file_name_, ',');
          size_t i = 0;
          for (const auto& mod : utils::split(out_mod_name_, ',')) {
            ParametersList outm;
            if (out_files.size() > i && !out_files.at(i).empty())
              outm.set<std::string>("filename", out_files.at(i));
            runParameters()->addEventExporter(EventExporterFactory::get().build(mod, outm));
            ++i;
          }
        }
        return *this;
      }
      /// Store a configuration into a LPAIR steering card
      void write(const std::string& filename) const override {
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

        std::ofstream f(filename, std::fstream::out | std::fstream::trunc);
        if (!f.is_open())
          throw CG_ERROR("LpairHandler") << "Failed to open file '" << filename << "' for writing.";
        for (const auto& ln : out_map)
          f << ln.second;
        f.close();
      }

    private:
      LpairHandler& setRunParameters(const RunParameters*) override;
      /// Single parameter handler
      /// \tparam T Parameter type
      template <typename T>
      struct Parameter {
        std::string key, description;
        T* value{nullptr};
      };
      /// Register a parameter to be steered to a configuration variable
      template <typename T>
      void registerParameter(const std::string& /*key*/, const std::string& /*description*/, T* /*def*/) {}
      template <typename T>
      void registerProcessParameter(const std::string& key,
                                    const std::string& description,
                                    const std::string& proc_key) {
        registerParameter<T>(key, description, &proc_params_->operator[]<T>(proc_key));
      }
      /// Register a kinematics block parameter to be steered
      template <typename T>
      void registerKinematicsParameter(const std::string& key,
                                       const std::string& description,
                                       const std::string& kin_key) {
        registerParameter<T>(key, description, &kin_params_->operator[]<T>(kin_key));
      }
      template <typename T>
      void registerGenerationParameter(const std::string& key,
                                       const std::string& description,
                                       const std::string& gen_key) {
        registerParameter<T>(key, description, &gen_params_->operator[]<T>(gen_key));
      }
      template <typename T>
      void registerIntegratorParameter(const std::string& key,
                                       const std::string& description,
                                       const std::string& int_key) {
        registerParameter<T>(key, description, &int_params_->operator[]<T>(int_key));
      }
      /// Set a parameter value
      template <typename T>
      inline void set(const std::string& /*key*/, const T& /*value*/) {}
      /// Retrieve a parameter value
      template <typename T>
      inline T get(const std::string& /*key*/) const {
        return T();
      }

      void setParameter(const std::string& key, const std::string& value);
      std::string parameter(std::string key) const;
      inline std::string describe(std::string key) const {
        if (p_strings_.count(key))
          return p_strings_.find(key)->second.description;
        if (p_ints_.count(key))
          return p_ints_.find(key)->second.description;
        if (p_doubles_.count(key))
          return p_doubles_.find(key)->second.description;
        return kInvalidStr;
      }

      static constexpr int kInvalidInt = -999999;
      static constexpr double kInvalidDbl = 999.999;
      static constexpr const char* kInvalidStr = "(null)";

      std::unordered_map<std::string, Parameter<std::string> > p_strings_;
      std::unordered_map<std::string, Parameter<double> > p_doubles_;
      std::unordered_map<std::string, Parameter<int> > p_ints_;

      void init();
      std::shared_ptr<ParametersList> proc_params_, kin_params_, gen_params_, int_params_;
      int timer_{0}, iend_{1}, log_level_{0}, ext_log_{0};
      int str_fun_{11}, sr_type_{1}, lepton_id_{0};
      std::string proc_name_, evt_mod_name_, out_mod_name_;
      std::string out_file_name_, addons_list_;
      std::string kmr_grid_path_, mstw_grid_path_, pdg_input_path_;
    };

    //----- specialised registerers

    /// Register a string parameter
    template <>
    inline void LpairHandler::registerParameter<std::string>(const std::string& key,
                                                             const std::string& description,
                                                             std::string* def) {
      p_strings_[key] = Parameter<std::string>{key, description, def};
    }
    /// Register a double floating point parameter
    template <>
    inline void LpairHandler::registerParameter<double>(const std::string& key,
                                                        const std::string& description,
                                                        double* def) {
      p_doubles_[key] = Parameter<double>{key, description, def};
    }
    /// Register an integer parameter
    template <>
    inline void LpairHandler::registerParameter<int>(const std::string& key, const std::string& description, int* def) {
      p_ints_[key] = Parameter<int>{key, description, def};
    }

    //----- specialised setters

    template <>
    inline void LpairHandler::set<std::string>(const std::string& key, const std::string& value) {
      if (p_strings_.count(key))
        *p_strings_.at(key).value = value;
    }
    template <>
    inline void LpairHandler::set<double>(const std::string& key, const double& value) {
      if (p_doubles_.count(key))
        *p_doubles_.at(key).value = value;
    }
    template <>
    inline void LpairHandler::set<int>(const std::string& key, const int& value) {
      if (p_ints_.count(key))
        *p_ints_.at(key).value = value;
    }

    //----- specialised getters

    /// Retrieve a string parameter value
    template <>
    inline std::string LpairHandler::get(const std::string& key) const {
      if (p_strings_.count(key))
        return *p_strings_.at(key).value;
      return kInvalidStr;
    }
    /// Retrieve a floating point parameter value
    template <>
    inline double LpairHandler::get(const std::string& key) const {
      if (p_doubles_.count(key))
        return *p_doubles_.at(key).value;
      return -kInvalidDbl;
    }
    /// Retrieve an integer parameter value
    template <>
    inline int LpairHandler::get(const std::string& key) const {
      if (p_ints_.count(key))
        return *p_ints_.at(key).value;
      return -kInvalidInt;
    }

    std::string LpairHandler::parameter(std::string key) const {
      if (auto var = get<double>(key); var != -kInvalidDbl)
        return std::to_string(var);
      if (auto var = get<int>(key); var != -kInvalidInt)
        return std::to_string(var);
      return get<std::string>(key);
    }

    void LpairHandler::init() {
      //-------------------------------------------------------------------------------------------
      // Process/integration/hadronisation parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<std::string>("PROC", "Process name to simulate", &proc_name_);
      registerParameter<std::string>("HADR", "Hadronisation algorithm", &evt_mod_name_);
      registerParameter<std::string>("EVMD", "Events modification algorithms", &evt_mod_name_);
      registerParameter<std::string>("OUTP", "Output module", &out_mod_name_);
      registerParameter<std::string>("OUTF", "Output file name", &out_file_name_);
      registerParameter<std::string>("ADDN", "Additional libraries to load", &addons_list_);
      registerIntegratorParameter<std::string>("ITYP", "Integration algorithm", MODULE_NAME);
      registerIntegratorParameter<int>("NTRT", "Smoothen the integrand", "treat");
      registerIntegratorParameter<int>("NCVG", "Number of function calls", "numFunctionCalls");
      registerIntegratorParameter<int>("ITVG", "Number of integration iterations", "iterations");
      registerIntegratorParameter<int>("SEED", "Random generator seed", "seed");

      //-------------------------------------------------------------------------------------------
      // General parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<int>("TIMR", "Enable the time ticker", &timer_);
      registerParameter<int>("IEND", "Generation type", &iend_);
      registerParameter<int>("DEBG", "Debugging verbosity", &log_level_);
      registerParameter<int>("LOGE", "Extended logging", &ext_log_);
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
      registerKinematicsParameter<int>("MODE", "Subprocess' mode", "mode");
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

    LpairHandler& LpairHandler::setRunParameters(const RunParameters* params) {
      Handler::setRunParameters(params);
      str_fun_ = runParameters()->kinematics().incomingBeams().structureFunctions().name<int>();
      sr_type_ = runParameters()->kinematics().incomingBeams().structureFunctions().get<int>("sigmaRatio");
      //kmr_grid_path_ = kmr::GluonGrid::get().path();
      //mstw_grid_path_ =
      //pdg_input_path_ =
      iend_ = (int)runParameters()->generation().enabled();
      log_level_ = (int)utils::Logger::get().level();
      ext_log_ = utils::Logger::get().extended();
      proc_name_ = runParameters()->processName();
      *proc_params_ += runParameters()->process().parameters();
      if (proc_params_->has<ParticleProperties>("pair"))
        proc_params_->set<int>("pair", proc_params_->get<ParticleProperties>("pair").pdgid);
      if (proc_name_ == "pptoff" || proc_name_ == "pptoll" /* legacy */)
        lepton_id_ = (runParameters()->process().parameters().get<int>("pair") - 11) / 2. + 1;
      {
        std::vector<std::string> evt_mod;
        std::transform(runParameters()->eventModifiersSequence().begin(),
                       runParameters()->eventModifiersSequence().end(),
                       std::back_inserter(evt_mod),
                       [](const auto& mod) { return mod->name(); });
        evt_mod_name_ = utils::merge(evt_mod, ",");
      }
      {
        std::vector<std::string> out_mod, out_mod_file;
        for (const auto& out : runParameters()->eventExportersSequence()) {
          out_mod.emplace_back(out->name());
          out_mod_file.emplace_back(out->parameters().get<std::string>("filename"));
        }
        out_mod_name_ = utils::merge(out_mod, ",");
        out_file_name_ = utils::merge(out_mod_file, ",");
      }
      timer_ = (runParameters()->timeKeeper() != nullptr);

      *kin_params_ += runParameters()->kinematics().fullParameters();
      *int_params_ += runParameters()->integrator();
      *gen_params_ += runParameters()->generation().parameters();
      init();
      return *this;
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
  }  // namespace card
}  // namespace cepgen
using cepgen::card::LpairHandler;
REGISTER_CARD_HANDLER(".card", LpairHandler);
