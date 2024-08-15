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
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

using std::string;

static constexpr int kInvalidInt = -999999;
static constexpr double kInvalidDbl = 999.999;
static constexpr const char* kInvalidStr = "(null)";

#define REGISTER_LPAIR_CONTENT_TYPE             \
  __TYPE_ENUM(int, p_ints_, -kInvalidInt)       \
  __TYPE_ENUM(double, p_doubles_, -kInvalidDbl) \
  __TYPE_ENUM(std::string, p_strings_, kInvalidStr)

namespace cepgen::card {
  /// LPAIR-like steering cards parser and writer
  class LpairHandler final : public Handler {
  public:
    /// Read a LPAIR steering card
    explicit LpairHandler(const ParametersList& params) : Handler(params) {}

    static ParametersDescription description() {
      auto desc = Handler::description();
      desc.setDescription("LPAIR-like cards parser");
      return desc;
    }

    inline LpairHandler& parseFile(const std::string& filename) override {
      if (!utils::fileExists(filename))
        throw CG_FATAL("LpairHandler:parseFile") << "Unable to locate steering card \"" << filename << "\".";
      if (const auto file_content = utils::split(utils::readFile(filename), '\n'); !file_content.empty())
        return parseCommands(file_content);
      CG_WARNING("LpairHandler:parseFile") << "Empty steering card.";
      return *this;
    }

    inline LpairHandler& parseCommands(const std::vector<std::string>& commands) override {
      std::ostringstream os;
      init();
      for (auto line : commands) {
        if (line = utils::trim(line); line.empty() || line[0] == '#')  // skip comments
          continue;
        const auto fields = utils::split(line, ' ', true);  // parse all fields

        if (fields.size() < 2) {  // skip all invalid lines
          CG_WARNING("LpairHandler:parseCommands") << "Invalid command read: '" << line << "'.";
          continue;
        }
        const auto key = fields.at(0), value = fields.at(1);
        setParameter(key, value);
        if (const auto descr = describe(key); descr != kInvalidStr)
          os << utils::format("\n\t>> %-8s %-25s (%s)", key.data(), parameter(key).data(), descr.data());
      }

      CG_INFO("LpairHandler:parseCommands")
          << "LPAIR configuration successfully loaded! Now parsing the following parameters:" << os.str() << ".";
      parse();
      return *this;
    }

    /// Store a configuration into a LPAIR steering card
    inline void write(const std::string& filename) const override {
      std::ofstream file(filename, std::fstream::out | std::fstream::trunc);
      if (!file.is_open())
        throw CG_ERROR("LpairHandler") << "Failed to open file '" << filename << "' for writing.";
      for (const auto& it : p_strings_)
        if (it.second.value && !it.second.value->empty())
          file << utils::format(
              "%-8s %-20s ! %s\n", it.first.data(), it.second.value->data(), it.second.description.data());
      for (const auto& it : p_ints_)
        if (it.second.value && *it.second.value != kInvalidInt)
          file << utils::format("%-8s %-20d ! %s\n", it.first.data(), *it.second.value, it.second.description.data());
      for (const auto& it : p_doubles_)
        if (it.second.value && *it.second.value != Limits::INVALID)
          file << utils::format("%-8s %-20e ! %s\n", it.first.data(), *it.second.value, it.second.description.data());
      file.close();
    }

  private:
    LpairHandler& setRunParameters(const RunParameters*) override;

    void init();
    inline void parse() {
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
      utils::Logger::get().setLevel(static_cast<utils::Logger::Level>(log_level_));
      utils::Logger::get().setExtended(ext_log_);

      {  //--- parse the structure functions code
        auto& kin_params = proc_params_.operator[]<ParametersList>("kinematics");
        bool beam1_elastic{true}, beam2_elastic{true};
        if (pmod_ == "1")
          kin_params.set<int>("beam1id", PDG::electron);
        else {
          kin_params.set<int>("beam1id", PDG::proton);
          if (pmod_ != "2") {
            beam1_elastic = false;
            auto sf_params = StructureFunctionsFactory::get().describeParameters(pmod_).parameters();
            sf_params.set<ParametersList>("sigmaRatio",
                                          SigmaRatiosFactory::get().describeParameters(sr_type_).parameters());
            if (pmod_ == "205" /* MSTWgrid */ && !mstw_grid_path_.empty())
              sf_params.set<std::string>("gridPath", mstw_grid_path_);
            kin_params.set("structureFunctions", sf_params);
          }
        }
        if (emod_ == "1")
          kin_params.set<int>("beam2id", PDG::electron);
        else {
          kin_params.set<int>("beam2id", PDG::proton);
          if (emod_ != "2") {
            beam2_elastic = false;
            if (emod_ != pmod_)
              CG_WARNING("LpairHandler") << "Incoming particles' modes are inconsistent: PMOD=" << pmod_
                                         << ", EMOD=" << emod_ << ".";
          }
        }
        auto& mode = kin_params.operator[]<int>("mode");
        mode = static_cast<int>(
            beam1_elastic
                ? (beam2_elastic ? mode::Kinematics::ElasticElastic : mode::Kinematics::ElasticInelastic)
                : (beam2_elastic ? mode::Kinematics::InelasticElastic : mode::Kinematics::InelasticInelastic));
      }

      //--- parse the process name
      if (!proc_name_.empty() || !proc_params_.empty()) {
        if (!runParameters()->hasProcess() && proc_name_.empty())
          throw CG_FATAL("LpairHandler") << "Process name not specified!";
        if (runParameters()->hasProcess() && runParameters()->process().name() == proc_name_)
          proc_params_ = runParameters()->process().parameters() + proc_params_;
        if (proc_name_ == "pptoff" && lepton_id_ != 0)  // backward-compatibility for PPtoLL cards
          proc_params_.set<int>("pair", PDG::electron + (lepton_id_ - 1) * 2);
        runParameters()->setProcess(ProcessFactory::get().build(proc_name_, proc_params_));
      }

      if (!int_params_.name().empty())
        runParameters()->integrator() += int_params_;
      runParameters()->generation().setParameters(gen_params_);

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
    }

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
    void registerProcessParameter(const std::string& key, const std::string& description, const std::string& proc_key) {
      registerParameter<T>(key, description, &proc_params_.operator[]<T>(proc_key));
    }
    /// Register a kinematics block parameter to be steered
    template <typename T>
    void registerKinematicsParameter(const std::string& key,
                                     const std::string& description,
                                     const std::string& kin_key) {
      registerParameter<T>(
          key, description, &proc_params_.operator[]<ParametersList>("kinematics").operator[]<T>(kin_key));
    }
    template <typename T>
    void registerGenerationParameter(const std::string& key,
                                     const std::string& description,
                                     const std::string& gen_key) {
      registerParameter<T>(key, description, &gen_params_.operator[]<T>(gen_key));
    }
    template <typename T>
    void registerIntegratorParameter(const std::string& key,
                                     const std::string& description,
                                     const std::string& int_key) {
      registerParameter<T>(key, description, &int_params_.operator[]<T>(int_key));
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
    std::string parameter(const std::string& key) const;
    inline std::string describe(const std::string& key) const {
#define __TYPE_ENUM(type, map, default_val) \
  if (map.count(key))                       \
    return map.at(key).description;
      REGISTER_LPAIR_CONTENT_TYPE
#undef __TYPE_ENUM
      return kInvalidStr;
    }

    ParametersList proc_params_, gen_params_, int_params_;
    int timer_{0}, iend_{1}, log_level_{static_cast<int>(utils::Logger::get().level())}, ext_log_{0};
    std::string emod_{"2"}, pmod_{"2"};
    int sr_type_{1}, lepton_id_{0};
    std::string proc_name_, evt_mod_name_, out_mod_name_;
    std::string out_file_name_, addons_list_;
    std::string kmr_grid_path_, mstw_grid_path_, pdg_input_path_;

#define __TYPE_ENUM(type, map_name, default_val) std::map<std::string, Parameter<type> > map_name;
    REGISTER_LPAIR_CONTENT_TYPE
#undef __TYPE_ENUM
  };

  //----- specialised registerers
#define __TYPE_ENUM(type, map, default_val)                                        \
  template <>                                                                      \
  inline void LpairHandler::registerParameter<type>(                               \
      const std::string& key, const std::string& description, type* def) {         \
    map[key] = Parameter<type>{key, description, def};                             \
  }                                                                                \
  template <>                                                                      \
  inline void LpairHandler::set<type>(const std::string& key, const type& value) { \
    if (map.count(key))                                                            \
      *map.at(key).value = value;                                                  \
  }                                                                                \
  template <>                                                                      \
  inline type LpairHandler::get(const std::string& key) const {                    \
    if (map.count(key))                                                            \
      return *map.at(key).value;                                                   \
    return default_val;                                                            \
  }
  REGISTER_LPAIR_CONTENT_TYPE
#undef __TYPE_ENUM

  std::string LpairHandler::parameter(const std::string& key) const {
#define __TYPE_ENUM(type, map, default_val)          \
  if (auto var = get<type>(key); var != default_val) \
    return utils::toString(var);
    REGISTER_LPAIR_CONTENT_TYPE
#undef __TYPE_ENUM
    CG_ERROR("LpairHandler:parameter") << "Failed to retrieve a parameter with key '" << key << "'.";
    return kInvalidStr;
  }

  void LpairHandler::init() {
    //-------------------------------------------------------------------------------------------
    // Process/integration/hadronisation parameters
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
    registerProcessParameter<int>("METH", "Computation method (kT-factorisation)", "method");
    registerProcessParameter<int>("IPOL", "Polarisation states to consider", "polarisationStates");

    //-------------------------------------------------------------------------------------------
    // Process kinematics parameters
    registerParameter<std::string>("KMRG", "KMR grid interpolation path", &kmr_grid_path_);
    registerParameter<std::string>("MGRD", "MSTW grid interpolation path", &mstw_grid_path_);
    registerParameter<std::string>("PDGI", "Input file for PDG information", &pdg_input_path_);
    registerParameter<std::string>("PMOD", "Outgoing primary particles' mode", &pmod_);
    registerParameter<std::string>("EMOD", "Outgoing primary particles' mode", &emod_);
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
    const auto beam_mode = runParameters()->kinematics().incomingBeams().mode();
    pmod_ = (beam_mode == mode::Kinematics::InelasticElastic || beam_mode == mode::Kinematics::InelasticInelastic)
                ? runParameters()->kinematics().incomingBeams().structureFunctions().name()
                : (std::abs(runParameters()->kinematics().incomingBeams().positive().integerPdgId()) == PDG::electron
                       ? "1"
                       : "2");
    emod_ = (beam_mode == mode::Kinematics::ElasticInelastic || beam_mode == mode::Kinematics::InelasticInelastic)
                ? runParameters()->kinematics().incomingBeams().structureFunctions().name()
                : (std::abs(runParameters()->kinematics().incomingBeams().negative().integerPdgId()) == PDG::electron
                       ? "1"
                       : "2");
    sr_type_ = runParameters()->kinematics().incomingBeams().structureFunctions().get<int>("sigmaRatio");
    //kmr_grid_path_ = kmr::GluonGrid::get().path();
    //mstw_grid_path_ =
    //pdg_input_path_ =
    iend_ = static_cast<int>(runParameters()->generation().enabled());
    log_level_ = static_cast<int>(utils::Logger::get().level());
    ext_log_ = utils::Logger::get().extended();
    proc_name_ = runParameters()->processName();
    proc_params_ += runParameters()->process().parameters();
    if (proc_params_.has<ParticleProperties>("pair"))
      proc_params_.set<int>("pair", proc_params_.get<ParticleProperties>("pair").pdgid);
    if (proc_name_ == "pptoff" || proc_name_ == "pptoll" /* legacy */)
      lepton_id_ = (runParameters()->process().parameters().get<int>("pair") - PDG::electron) / 2. + 1;
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

    int_params_ += runParameters()->integrator();
    gen_params_ += runParameters()->generation().parameters();
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
}  // namespace cepgen::card
using cepgen::card::LpairHandler;
REGISTER_CARD_HANDLER(".card", LpairHandler);
