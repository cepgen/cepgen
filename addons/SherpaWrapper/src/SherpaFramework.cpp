/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

// CepGen includes
#include <CepGen/Core/Exception.h>
#include <CepGen/Utils/Message.h>
#include <CepGen/Utils/StreamCollector.h>
#include <CepGen/Utils/String.h>
// Sherpa includes
#include <ATOOLS/Math/Random.H>
#include <ATOOLS/Org/Exception.H>
#include <ATOOLS/Org/Library_Loader.H>
#include <ATOOLS/Org/Message.H>
#include <ATOOLS/Org/My_MPI.H>
#include <ATOOLS/Org/Run_Parameter.H>
#include <ATOOLS/Org/Scoped_Settings.H>
#include <ATOOLS/Org/Settings.H>
#include <ATOOLS/Org/Terminator_Objects.H>
#include <ATOOLS/Phys/Flavour.H>
#include <ATOOLS/Phys/Fragmentation_Base.H>
#include <PDF/Main/PDF_Base.H>
#include <SHERPA/Initialization/Initialization_Handler.H>
#include <SHERPA/Single_Events/Hadron_Decays.H>
#include <SHERPA/Single_Events/Hadronization.H>
#include <SHERPA/Single_Events/Jet_Evolution.H>
#include <unistd.h>

#include <array>
#include <csignal>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "CepGenSherpa/SherpaFramework.h"

namespace cepgen::sherpa {
  namespace {
    /// Preserve the signal and termination handlers of the host program, as the Sherpa exceptions handler overrides
    /// them at its construction
    class HostHandlersGuard {
    public:
      HostHandlersGuard() : terminate_(std::get_terminate()) {
        for (size_t i = 0; i < kSignals.size(); ++i)
          sigaction(kSignals.at(i), nullptr, &actions_.at(i));
      }
      ~HostHandlersGuard() {
        for (size_t i = 0; i < kSignals.size(); ++i)
          sigaction(kSignals.at(i), &actions_.at(i), nullptr);
        std::set_terminate(terminate_);
      }

    private:
      static constexpr std::array<int, 8> kSignals{SIGSEGV, SIGINT, SIGPIPE, SIGBUS, SIGFPE, SIGABRT, SIGTERM, SIGXCPU};
      std::array<struct sigaction, kSignals.size()> actions_{};
      const std::terminate_handler terminate_;
    };

    std::string describe(const ATOOLS::Exception& exception) {
      std::ostringstream os;
      os << exception;
      return os.str();
    }

    /// Temporary working directory for the Sherpa initialisation, holding the defaults run card, and receiving the
    /// (unused) process libraries directory structure; removed once the framework is set up
    class TemporaryDirectory {
    public:
      explicit TemporaryDirectory(const std::string& defaults)
          : path_(std::filesystem::temp_directory_path() / ("cepgen_sherpa_" + std::to_string(::getpid()))) {
        std::filesystem::create_directories(path_);
        std::ofstream(runCard()) << defaults;
      }
      ~TemporaryDirectory() {
        std::error_code ec;
        std::filesystem::remove_all(path_, ec);
      }
      std::string path() const { return path_.string(); }
      std::string runCard() const { return (path_ / "defaults.yaml").string(); }

    private:
      const std::filesystem::path path_;
    };

    /// Command-line arguments for the Sherpa settings handler
    class Arguments {
    public:
      void add(const std::string& argument) {
        storage_.emplace_back(argument);
        pointers_.clear();
      }
      int argc() const { return static_cast<int>(storage_.size()); }
      char** argv() {
        pointers_.clear();
        for (auto& argument : storage_)
          pointers_.emplace_back(argument.data());
        pointers_.emplace_back(nullptr);
        return pointers_.data();
      }

    private:
      std::vector<std::string> storage_{"cepgen"};
      std::vector<char*> pointers_;
    };

    /// Defaults otherwise registered by the SHERPA::Sherpa object, and required by the initialisation handler
    void registerRunDefaults() {
      auto& settings = ATOOLS::Settings::GetMainSettings();
      settings["PRINT_VERSION_INFO"].SetDefault(false);
      settings["INIT_ONLY"].SetDefault(0);
      settings["SHOW_NTRIALS"].SetDefault(false);
      settings["DEBUG_INTERVAL"].SetDefault(0);
      settings["DEBUG_STEP"].SetDefault(-1);
      settings["EVENT_DISPLAY_INTERVAL"].SetDefault(100);
      settings["EVT_OUTPUT"].SetDefault(ATOOLS::msg->Level());
      settings["MSG_LIMIT"].SetDefault(20);
      ATOOLS::msg->SetLimit(settings["MSG_LIMIT"].Get<int>());
      settings["EVT_OUTPUT_START"].SetDefault(0);
    }

    void addRunCard(Arguments& arguments, const std::string& path) {
      // file names are recognised by the Sherpa command line parser by the absence of any ':' or '=' character
      if (path.find(':') != std::string::npos || path.find('=') != std::string::npos)
        throw CG_FATAL("sherpa:Framework")
            << "Sherpa run card path '" << path << "' cannot contain any ':' or '=' character.";
      arguments.add(std::filesystem::absolute(path).string());
    }
  }  // namespace

  bool Framework::Configuration::operator==(const Configuration& oth) const {
    return arguments == oth.arguments && run_card == oth.run_card && defaults == oth.defaults;
  }

  Framework& Framework::get() {  // standard Meyer's singleton
    static Framework framework;
    return framework;
  }

  void Framework::initialise(const Configuration& configuration) {
    std::lock_guard<std::mutex> lock(init_mutex_);
    if (initialised_) {
      if (configuration != configuration_)
        throw CG_FATAL("sherpa:Framework")
            << "Sherpa was already set up in this process with a different configuration.\n\t"
            << "Current: " << utils::merge(configuration_.arguments, ", ") << " (run card: '" << configuration_.run_card
            << "')\n\t" << "Requested: " << utils::merge(configuration.arguments, ", ") << " (run card: '"
            << configuration.run_card << "').";
      return;
    }
    if (!configuration.run_card.empty() && !std::filesystem::exists(configuration.run_card))
      throw CG_FATAL("sherpa:Framework") << "Sherpa run card '" << configuration.run_card << "' does not exist.";

    const TemporaryDirectory working_directory(configuration.defaults);
    Arguments arguments;
    arguments.add("SHERPA_CPP_PATH: " + working_directory.path());  // do not pollute the current directory
    for (const auto& argument : configuration.arguments)
      arguments.add(argument);
    addRunCard(arguments, working_directory.runCard());  // run cards are listed by increasing priority
    if (!configuration.run_card.empty())
      addRunCard(arguments, configuration.run_card);

    std::string output, error;
    bool success = false;
    try {
      const HostHandlersGuard guard;
      utils::StreamCollector collector(output);
      ATOOLS::mpi = new ATOOLS::My_MPI();
      ATOOLS::exh = new ATOOLS::Terminator_Object_Handler();
      ATOOLS::msg = new ATOOLS::Message();
      ATOOLS::rpa = new ATOOLS::Run_Parameter();
      ATOOLS::Settings::InitializeMainSettings(arguments.argc(), arguments.argv());
      ATOOLS::ran = new ATOOLS::Random(1234);  // overwritten
      ATOOLS::s_loader = new ATOOLS::Library_Loader();
      PDF::pdfdefs = new PDF::PDF_Defaults();
      init_handler_ = new SHERPA::Initialization_Handler();
      registerRunDefaults();
      success = init_handler_->InitializeTheFramework();
    } catch (const ATOOLS::Exception& exception) {
      error = describe(exception);
    } catch (const std::exception& exception) {
      error = exception.what();
    }
    if (!success)
      throw CG_FATAL("sherpa:Framework") << "Failed to initialise the Sherpa framework"
                                         << (error.empty() ? "" : ": " + error) << ".\nSherpa output:\n"
                                         << output;
    CG_DEBUG("sherpa:Framework") << "Sherpa initialisation output:\n" << output;

    fragmentation_ = init_handler_->GetFragmentation();
    if (!fragmentation_ || fragmentation_->Name() == "None")
      throw CG_FATAL("sherpa:Framework") << "No fragmentation algorithm was set up by Sherpa.";
    hadron_decays_ = init_handler_->GetHDHandler() != nullptr;
    secondary_phases_ = {new SHERPA::Jet_Evolution(init_handler_->GetMatrixElementHandler(),
                                                   init_handler_->GetHardDecayHandler(),
                                                   init_handler_->GetHDHandler(),
                                                   init_handler_->GetMIHandlers(),
                                                   init_handler_->GetSoftCollisionHandlers(),
                                                   init_handler_->GetShowerHandlers(),
                                                   init_handler_->GetRemnantHandlers()),
                         new SHERPA::Hadronization(init_handler_->GetColourReconnectionHandler(), fragmentation_),
                         new SHERPA::Hadron_Decays(init_handler_->GetHDHandler())};
    CG_INFO("sherpa:Framework") << "Sherpa framework initialised with '" << fragmentation_->Name()
                                << "' fragmentation algorithm, " << (hadron_decays_ ? "with" : "without")
                                << " hadron decays.";
    configuration_ = configuration;
    initialised_ = true;
  }

  ATOOLS::Fragmentation_Base& Framework::fragmentation() const {
    if (!fragmentation_)
      throw CG_FATAL("sherpa:Framework") << "Sherpa framework is not initialised.";
    return *fragmentation_;
  }

  bool Framework::hadronDecaysEnabled() const { return hadron_decays_; }

  double Framework::partonMass(long pdg_id, bool constituent) {
    const ATOOLS::Flavour flavour(static_cast<kf_code>(std::labs(pdg_id)), pdg_id < 0);
    return constituent ? flavour.HadMass() : flavour.Mass(true);
  }
}  // namespace cepgen::sherpa
