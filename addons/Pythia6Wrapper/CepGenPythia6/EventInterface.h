/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2026  Laurent Forthomme
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

#ifndef CepGenPythia6_EventInterface_h
#define CepGenPythia6_EventInterface_h

#include <vector>

#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/Hadroniser.h"

namespace cepgen {
  class Event;
}  // namespace cepgen

namespace cepgen::pythia6 {
  /// Interface between the CepGen and the Pythia 6 event contents
  class EventInterface {
  public:
    /// Build an interface to a CepGen event
    /// \param[in] event Event to be modified (not owned)
    /// \note Dissociated beams are left untouched by decay(), and are to be fragmented through a DiffractiveState
    explicit EventInterface(Event& event);

    bool decay();  ///< Run Pythia 6 decay/fragmentation algorithm on full event content
    size_t numStrings() const { return evt_strings_.size(); }  ///< Number of string-bound partons in the event

    /// Attach the hadrons from a dissociated beam fragmentation to its parent system
    /// \param[in] role Role of the dissociated system in the event
    /// \param[in] fragments Hadrons produced, in the rest frame of the dissociated system
    void update(Particle::Role role, const hadr::Hadroniser::Fragments& fragments);

  private:
    void fillEventBlock();

    Event& cepgen_event_;  ///< CepGen event content (not owning)

    using string_t = std::vector<int>;
    std::vector<string_t> evt_strings_;
  };
}  // namespace cepgen::pythia6

#endif
