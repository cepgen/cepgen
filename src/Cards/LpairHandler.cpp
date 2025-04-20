/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

using namespace cepgen;
using namespace cepgen::card;
using namespace std::string_literals;
using std::string;

static constexpr int kInvalidInt = -999999;
static constexpr double kInvalidDbl = 999.999;
static constexpr const char* kInvalidStr = "(null)";

#define REGISTER_LPAIR_CONTENT_TYPE             \
  __TYPE_ENUM(int, p_ints_, -kInvalidInt)       \
  __TYPE_ENUM(double, p_doubles_, -kInvalidDbl) \
  __TYPE_ENUM(std::string, p_strings_, kInvalidStr)

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

  LpairHandler& parseFile(const std::string& filename) override {
    if (!utils::fileExists(filename))
      throw CG_FATAL("LpairHandler:parseFile") << "Unable to locate steering card \"" << filename << "\".";
    if (const auto file_content = utils::split(utils::readFile(filename), '\n'); !file_content.empty())
      return parseCommands(file_content);
    CG_WARNING("LpairHandler:parseFile") << "Empty steering card.";
    return *this;
  }

  LpairHandler& parseCommands(const std::vector<std::string>& commands) override {
    std::ostringstream os;
    init();
    for (auto line : commands) {
      if (line = utils::trim(line); line.empty() || line[0] == '#')  // skip comments
        continue;
      const auto fields = utils::split(line, ' ', true);  // parse all fields
      if (fields.size() < 2) {                            // skip all invalid lines
        CG_WARNING("LpairHandler:parseCommands") << "Invalid command read: '" << line << "'.";
        continue;
      }
      const auto &key = fields.at(0), &value = fields.at(1);
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
  void write(const std::string& filename) const override {
    std::ofstream file(filename, std::fstream::out | std::fstream::trunc);
    if (!file.is_open())
      throw CG_ERROR("LpairHandler") << "Failed to open file '" << filename << "' for writing.";
    for (const auto& [key, value] : p_strings_)
      if (value.value && !value.value->empty())
        file << utils::format("%-8s %-20s ! %s\n", key.data(), value.value->data(), value.description.data());
    for (const auto& [key, value] : p_ints_)
      if (value.value && *value.value != kInvalidInt)
        file << utils::format("%-8s %-20d ! %s\n", key.data(), *value.value, value.description.data());
    for (const auto& [key, value] : p_doubles_)
      if (value.value && *value.value != Limits::INVALID)
        file << utils::format("%-8s %-20e ! %s\n", key.data(), *value.value, value.description.data());
    file.close();
  }

private:
  LpairHandler& setRunParameters(const RunParameters*) override;

  void init();
  void parse() {
    for (const auto& lib : utils::split(addons_list_, ','))
      loadLibrary(lib);

    // parse the PDG library
    if (!pdg_input_path_.empty())
      pdg::MCDFileParser::parse(pdg_input_path_);
    if (!kmr_grid_path_.empty())
      kmr::GluonGrid::get(ParametersList().set("path", kmr_grid_path_));

    if (timer_)  // build the ticker if required
      runParameters()->setTimeKeeper(new utils::TimeKeeper);
    utils::Logger::get().setLevel(static_cast<utils::Logger::Level>(log_level_));
    utils::Logger::get().setExtended(ext_log_);

    {  // parse the structure functions code
      auto& kinematics_params = proc_params_.operator[]<ParametersList>("kinematics");
      bool beam1_elastic{true}, beam2_elastic{true};
      if (p_mode_ == "1")
        kinematics_params.set<int>("beam1id", PDG::electron);
      else {
        kinematics_params.set<int>("beam1id", PDG::proton);
        if (p_mode_ != "2") {
          beam1_elastic = false;
          auto sf_params = StructureFunctionsFactory::get().describeParameters(p_mode_).parameters();
          sf_params.set("sigmaRatio", SigmaRatiosFactory::get().describeParameters(sr_type_).parameters());
          if (p_mode_ == "205" /* MSTW grid */ && !mstw_grid_path_.empty())
            sf_params.set("gridPath", mstw_grid_path_);
          kinematics_params.set("structureFunctions", sf_params);
        }
      }
      if (e_mode_ == "1")
        kinematics_params.set<int>("beam2id", PDG::electron);
      else {
        kinematics_params.set<int>("beam2id", PDG::proton);
        if (e_mode_ != "2") {
          beam2_elastic = false;
          if (e_mode_ != p_mode_)
            CG_WARNING("LpairHandler") << "Incoming particles' modes are inconsistent: PMOD=" << p_mode_ << ", "
                                       << "EMOD"s << "=" << e_mode_ << ".";
        }
      }
      auto& mode = kinematics_params.operator[]<int>("mode");
      mode = static_cast<int>(
          beam1_elastic ? (beam2_elastic ? mode::Kinematics::ElasticElastic : mode::Kinematics::ElasticInelastic)
                        : (beam2_elastic ? mode::Kinematics::InelasticElastic : mode::Kinematics::InelasticInelastic));
    }

    if (!proc_name_.empty() || !proc_params_.empty()) {  // parse the process parameters
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

    for (const auto& mod : utils::split(evt_mod_name_, ','))  // parse the hadronisation algorithm name
      runParameters()->addModifier(EventModifierFactory::get().build(mod, ParametersList()));

    if (!out_mod_name_.empty()) {  // parse the output module name
      const auto& out_files = utils::split(out_file_name_, ',');
      size_t i = 0;
      for (const auto& mod : utils::split(out_mod_name_, ',')) {
        ParametersList output_params;
        if (out_files.size() > i && !out_files.at(i).empty())
          output_params.set("filename", out_files.at(i));
        runParameters()->addEventExporter(EventExporterFactory::get().build(mod, output_params));
        ++i;
      }
    }
  }

  /// Single parameter handler
  /// \tparam T Parameter type
  template <typename T>
  struct Parameter {
    std::string key;
    std::string description;
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
  void registerKinematicsParameter(const std::string& key, const std::string& description, const std::string& kin_key) {
    registerParameter<T>(
        key, description, &proc_params_.operator[]<ParametersList>("kinematics").operator[]<T>(kin_key));
  }
  template <typename T>
  void registerGenerationParameter(const std::string& key, const std::string& description, const std::string& gen_key) {
    registerParameter<T>(key, description, &gen_params_.operator[]<T>(gen_key));
  }
  template <typename T>
  void registerIntegratorParameter(const std::string& key, const std::string& description, const std::string& int_key) {
    registerParameter<T>(key, description, &int_params_.operator[]<T>(int_key));
  }
  /// Set a parameter value
  template <typename T>
  void set(const std::string& /*key*/, const T& /*value*/) {}
  /// Retrieve a parameter value
  template <typename T>
  T get(const std::string& /*key*/) const {
    return T();
  }

  void setParameter(const std::string& key, const std::string& value);
  std::string parameter(const std::string& key) const;
  std::string describe(const std::string& key) const {
#define __TYPE_ENUM(type, map, default_val) \
  if (map.count(key))                       \
    return map.at(key).description;
    REGISTER_LPAIR_CONTENT_TYPE
#undef __TYPE_ENUM
    return kInvalidStr;
  }

  ParametersList proc_params_, gen_params_, int_params_;
  int timer_{0}, generation_mode_{1}, log_level_{static_cast<int>(utils::Logger::get().level())}, ext_log_{0};
  std::string e_mode_{"2"}, p_mode_{"2"};
  int sr_type_{1}, lepton_id_{0};
  std::string proc_name_, evt_mod_name_, out_mod_name_;
  std::string out_file_name_, addons_list_;
  std::string kmr_grid_path_, mstw_grid_path_, pdg_input_path_;

#define __TYPE_ENUM(type, map_name, default_val) std::map<std::string, Parameter<type> > map_name;
  REGISTER_LPAIR_CONTENT_TYPE
#undef __TYPE_ENUM
};

// specialised registers
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
  registerParameter<std::string>("PROC"s, "Process name to simulate", &proc_name_);
  registerParameter<std::string>("HADR"s, "Hadronisation algorithm", &evt_mod_name_);
  registerParameter<std::string>("EVMD"s, "Events modification algorithms", &evt_mod_name_);
  registerParameter<std::string>("OUTP"s, "Output module", &out_mod_name_);
  registerParameter<std::string>("OUTF"s, "Output file name", &out_file_name_);
  registerParameter<std::string>("ADDN"s, "Additional libraries to load", &addons_list_);
  registerIntegratorParameter<std::string>("ITYP"s, "Integration algorithm", MODULE_NAME);
  registerIntegratorParameter<int>("NTRT"s, "Smoothen the integrand", "treat");
  registerIntegratorParameter<int>("NCVG"s, "Number of function calls", "numFunctionCalls");
  registerIntegratorParameter<int>("ITVG"s, "Number of integration iterations", "iterations");
  registerIntegratorParameter<int>("SEED"s, "Random generator seed", "seed");

  //-------------------------------------------------------------------------------------------
  // General parameters
  registerParameter("TIMR"s, "Enable the time ticker", &timer_);
  registerParameter("IEND"s, "Generation type", &generation_mode_);
  registerParameter("DEBG"s, "Debugging verbosity", &log_level_);
  registerParameter("LOGE"s, "Extended logging", &ext_log_);
  registerGenerationParameter<int>("NTHR"s, "Number of threads to use for events generation", "numThreads"s);
  registerGenerationParameter<int>("NCSG"s, "Number of points to probe", "numPoints"s);
  registerGenerationParameter<int>("NGEN"s, "Number of events to generate", "maxgen"s);
  registerGenerationParameter<int>("NPRN"s, "Number of events before printout", "printEvery"s);

  //-------------------------------------------------------------------------------------------
  // Process-specific parameters
  registerProcessParameter<int>("METH"s, "Computation method (kT-factorisation)", "method"s);
  registerProcessParameter<int>("IPOL"s, "Polarisation states to consider", "polarisationStates"s);

  //-------------------------------------------------------------------------------------------
  // Process kinematics parameters
  registerParameter("KMRG"s, "KMR grid interpolation path", &kmr_grid_path_);
  registerParameter("MGRD"s, "MSTW grid interpolation path", &mstw_grid_path_);
  registerParameter("PDGI"s, "Input file for PDG information", &pdg_input_path_);
  registerParameter("PMOD"s, "Outgoing primary particles' mode", &p_mode_);
  registerParameter("EMOD"s, "Outgoing primary particles' mode", &e_mode_);
  registerParameter("RTYP"s, "R-ratio computation type", &sr_type_);
  registerProcessParameter<int>("PAIR"s, "Outgoing particles' PDG id", "pair"s);
  registerKinematicsParameter<std::string>("FFAC"s, "Form factors for the incoming beams", "formFactors"s);
  registerKinematicsParameter<int>("MODE"s, "Subprocess' mode", "mode"s);
  registerKinematicsParameter<int>("INA1"s, "Heavy ion atomic weight (1st incoming beam)", "beam1A"s);
  registerKinematicsParameter<int>("INZ1"s, "Heavy ion atomic number (1st incoming beam)", "beam1Z"s);
  registerKinematicsParameter<int>("INA2"s, "Heavy ion atomic weight (2nd incoming beam)", "beam2A"s);
  registerKinematicsParameter<int>("INZ2"s, "Heavy ion atomic number (2nd incoming beam)", "beam2Z"s);
  registerKinematicsParameter<double>("INP1"s, "Momentum (1st primary particle)", "beam1pz"s);
  registerKinematicsParameter<double>("INP2"s, "Momentum (2nd primary particle)", "beam2pz"s);
  registerKinematicsParameter<double>("INPP"s, "Momentum (1st primary particle)", "beam1pz"s);
  registerKinematicsParameter<double>("INPE"s, "Momentum (2nd primary particle)", "beam2pz"s);
  registerKinematicsParameter<double>(
      "PTCT"s, "Minimal transverse momentum (single central outgoing particle)", "ptmin"s);
  registerKinematicsParameter<double>(
      "PTMX"s, "Maximal transverse momentum (single central outgoing particle)", "ptmax"s);
  registerKinematicsParameter<double>("MSCT"s, "Minimal central system mass", "invmassmin"s);
  registerKinematicsParameter<double>("MSMX"s, "Maximal central system mass", "invmassmax"s);
  registerKinematicsParameter<double>("ECUT"s, "Minimal energy (single central outgoing particle)", "energysummin"s);
  registerKinematicsParameter<double>("ETMN"s, "Minimal pseudo-rapidity (central outgoing particles)", "etamin"s);
  registerKinematicsParameter<double>("ETMX"s, "Maximal pseudo-rapidity (central outgoing particles)", "etamax"s);
  registerKinematicsParameter<double>("YMIN"s, "Minimal rapidity (central outgoing particles)", "rapiditymin"s);
  registerKinematicsParameter<double>("YMAX"s, "Maximal rapidity (central outgoing particles)", "rapiditymax"s);
  registerKinematicsParameter<double>(
      "PDMN"s, "Minimal transverse momentum difference (central outgoing particles)", "ptdiffmin"s);
  registerKinematicsParameter<double>(
      "PDMX"s, "Maximal transverse momentum difference (central outgoing particles)", "ptdiffmax"s);
  registerKinematicsParameter<double>("Q2MN"s, "Minimal Q^2 = -q^2 (exchanged parton)", "q2min"s);
  registerKinematicsParameter<double>("Q2MX"s, "Maximal Q^2 = -q^2 (exchanged parton)", "q2max"s);
  registerKinematicsParameter<double>("QTMN"s, "Minimal Q_T (exchanged parton)", "qtmin"s);
  registerKinematicsParameter<double>("QTMX"s, "Maximal Q_T (exchanged parton)", "qtmax"s);
  registerKinematicsParameter<double>("MXMN"s, "Minimal invariant mass of proton remnants", "mxmin"s);
  registerKinematicsParameter<double>("MXMX"s, "Maximal invariant mass of proton remnants", "mxmax"s);
  registerKinematicsParameter<double>("XIMN"s, "Minimal fractional momentum loss of outgoing proton (xi)", "ximin"s);
  registerKinematicsParameter<double>("XIMX"s, "Maximal fractional momentum loss of outgoing proton (xi)", "ximax"s);
  registerKinematicsParameter<double>("YJMN"s, "Minimal remnant jet rapidity", "yjmin"s);
  registerKinematicsParameter<double>("YJMX"s, "Maximal remnant jet rapidity", "yjmax"s);

  //-------------------------------------------------------------------------------------------
  // PPtoLL cards backward compatibility
  registerIntegratorParameter<int>("NTREAT"s, "Smoothen the integrand", "treat"s);
  registerIntegratorParameter<int>("ITMX"s, "Number of integration iterations", "iterations"s);
  registerIntegratorParameter<int>("NCVG"s, "Number of function calls to perform", "numFunctionCalls"s);
  registerProcessParameter<int>("METHOD"s, "Computation method (kT-factorisation)", "method"s);
  registerParameter("LEPTON"s, "Outgoing leptons' flavour", &lepton_id_);
  registerKinematicsParameter<double>(
      "PTMIN"s, "Minimal transverse momentum (single central outgoing particle)", "ptmin"s);
  registerKinematicsParameter<double>(
      "PTMAX"s, "Maximal transverse momentum (single central outgoing particle)", "ptmax"s);
  registerKinematicsParameter<double>("Q1TMIN"s, "Minimal Q_T (exchanged parton)", "qtmin"s);
  registerKinematicsParameter<double>("Q1TMAX"s, "Maximal Q_T (exchanged parton)", "qtmax"s);
  registerKinematicsParameter<double>("Q2TMIN"s, "Minimal Q_T (exchanged parton)", "qtmin"s);
  registerKinematicsParameter<double>("Q2TMAX"s, "Maximal Q_T (exchanged parton)", "qtmax"s);
  registerKinematicsParameter<double>("MXMIN"s, "Minimal invariant mass of proton remnants", "mxmin"s);
  registerKinematicsParameter<double>("MXMAX"s, "Maximal invariant mass of proton remnants", "mxmax"s);
}

LpairHandler& LpairHandler::setRunParameters(const RunParameters* params) {
  Handler::setRunParameters(params);
  const auto beam_mode = runParameters()->kinematics().incomingBeams().mode();
  p_mode_ =
      (beam_mode == mode::Kinematics::InelasticElastic || beam_mode == mode::Kinematics::InelasticInelastic)
          ? runParameters()->kinematics().incomingBeams().structureFunctions().name()
          : (std::abs(runParameters()->kinematics().incomingBeams().positive().integerPdgId()) == PDG::electron ? "1"
                                                                                                                : "2");
  e_mode_ =
      (beam_mode == mode::Kinematics::ElasticInelastic || beam_mode == mode::Kinematics::InelasticInelastic)
          ? runParameters()->kinematics().incomingBeams().structureFunctions().name()
          : (std::abs(runParameters()->kinematics().incomingBeams().negative().integerPdgId()) == PDG::electron ? "1"
                                                                                                                : "2");
  sr_type_ = runParameters()->kinematics().incomingBeams().structureFunctions().get<int>("sigmaRatio");
  //kmr_grid_path_ = kmr::GluonGrid::get().path();
  //mstw_grid_path_ =
  //pdg_input_path_ =
  generation_mode_ = static_cast<int>(runParameters()->generation().enabled());
  log_level_ = static_cast<int>(utils::Logger::get().level());
  ext_log_ = utils::Logger::get().extended();
  proc_name_ = runParameters()->processName();
  proc_params_ += runParameters()->process().parameters();
  if (proc_params_.has<ParticleProperties>("pair"))
    proc_params_.set<int>("pair", proc_params_.get<ParticleProperties>("pair").pdgid);
  if (proc_name_ == "pptoff"s || proc_name_ == "pptoll"s /* legacy */)
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
      set(key, std::stod(value));
      return;
    } catch (const std::logic_error&) {
      for (const auto& let : value)
        if (isalpha(let) && let != 'E' && let != 'e') {
          set(key, value);
          return;
        }
      throw CG_FATAL("LpairHandler:setParameter")
          << "Failed to parse a floating-point parameter \"" << key << "\" → \"" << value << "\"!";
    }
  try {
    set(key, std::stoi(value));
  } catch (const std::logic_error&) {
    try {
      set(key, value);
    } catch (const std::logic_error&) {
      throw CG_FATAL("LpairHandler:setParameter")
          << "Failed to add the parameter \"" << key << "\" → \"" << value << "\"!";
    }
  }
}
REGISTER_CARD_HANDLER(".card", LpairHandler);
