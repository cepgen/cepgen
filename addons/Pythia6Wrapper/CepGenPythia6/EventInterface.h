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

#ifndef CepGenPythia6_EventInterface_h
#define CepGenPythia6_EventInterface_h

#include <vector>

#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/Modes.h"

namespace cepgen {
  class Event;
  namespace utils {
    class RandomGenerator;
  }
}  // namespace cepgen

namespace cepgen::pythia6 {
  /// Interface to the Pythia 6 event content.
  class EventInterface {
  public:
    explicit EventInterface(Event&, mode::Kinematics, utils::RandomGenerator*);

    void prepareHadronisation();  ///< Add/edit event attributes to prepare for its fragmentation/hadronisation
    size_t numStrings() const { return evt_strings_.size(); }  ///< Number of string-bound partons in the event
    void run();                                                ///< Run the fragmentation/hadronisation algorithm

  private:
    void fillEventBlock();

    Event& evt_;                   // NOT owning
    utils::RandomGenerator* rnd_;  ///< Random number generator engine (not owning)
    std::vector<Particle::Role> roles_;

    std::pair<short, short> pickPartonsContent() const;

    using string_t = std::vector<int>;
    std::vector<string_t> evt_strings_;
  };
}  // namespace cepgen::pythia6

#endif
