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

#ifndef CepGenSherpa_SherpaFramework_h
#define CepGenSherpa_SherpaFramework_h

#include <mutex>
#include <string>
#include <vector>

namespace ATOOLS {
  class Fragmentation_Base;
}  // namespace ATOOLS
namespace SHERPA {
  class Initialization_Handler;
  class Event_Phase_Handler;
}  // namespace SHERPA

namespace cepgen::sherpa {
  /// Stripped Sherpa framework
  class Framework {
  public:
    /// Sherpa configuration
    struct Configuration {
      std::vector<std::string> arguments;
      std::string run_card;  ///< Optional user YAML run card
      std::string defaults;  ///< YAML block of defaults
      bool operator==(const Configuration&) const;
      bool operator!=(const Configuration& oth) const { return !operator==(oth); }
    };

    static Framework& get();  ///< Retrieve the process-wide instance

    /// Boot the Sherpa framework
    /// \note A second call with an identical configuration is a no-op, a different one is a fatal error
    void initialise(const Configuration&);
    bool initialised() const { return initialised_; }

    ATOOLS::Fragmentation_Base& fragmentation() const;  ///< Fragmentation algorithm
    bool hadronDecaysEnabled() const;                   ///< Decay unstable hadron?
    /// Event phases treating fragmentation products (hadron decays products showering, partons hadronisation, hadrons decays)
    const std::vector<SHERPA::Event_Phase_Handler*>& secondaryPhases() const { return secondary_phases_; }
    std::mutex& mutex() { return mutex_; }

    /// Mass of a parton, as defined in the Sherpa particles table
    /// \param[in] pdg_id Signed PDG identifier of the parton
    /// \param[in] constituent Constituent (hadronisation) mass, instead of the perturbative one
    static double partonMass(long pdg_id, bool constituent);

  private:
    Framework() = default;  // make this a singleton

    Configuration configuration_;
    SHERPA::Initialization_Handler* init_handler_{nullptr};
    std::vector<SHERPA::Event_Phase_Handler*> secondary_phases_;
    ATOOLS::Fragmentation_Base* fragmentation_{nullptr};  //NOT owning
    bool hadron_decays_{false};
    std::mutex init_mutex_;
    std::mutex mutex_;
    bool initialised_{false};
  };
}  // namespace cepgen::sherpa

#endif
