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

#ifndef CepGenPythia8_EventInterface_h
#define CepGenPythia8_EventInterface_h

#include <CepGen/Event/Particle.h>
#include <CepGen/Physics/DiffractiveState.h>

#include <unordered_map>

namespace Pythia8 {
  class Particle;
  class Pythia;
}  // namespace Pythia8
namespace cepgen {
  class Event;
}  // namespace cepgen

namespace cepgen::pythia8 {
  /// Interface between CepGen/Pythia event contents
  class EventInterface {
  public:
    /// Build an interface to a CepGen event
    /// \param[in] event Event to be modified (NOT owned)
    /// \param[in] pythia Pythia 8 engine, initialised with ProcessLevel:all = off (NOT owned)
    /// \note Dissociated beams are left untouched by decay(), and are to be fragmented through a DiffractiveState
    explicit EventInterface(Event& event, Pythia8::Pythia& pythia);

    /// Run the Pythia 8 decay/fragmentation algorithm on the central system
    /// \param[out] weight Product of the open decay channels fractions for all particles decayed
    bool decay(double& weight);

    /// Attach the hadrons from a dissociated beam fragmentation to its parent system
    /// \param[in] role Role of the dissociated system in the event
    /// \param[in] fragments Hadrons produced, in the rest frame of the dissociated system
    void update(Particle::Role role, const hadr::Fragments& fragments);

    static void registerParticle(const Pythia8::Particle&);  ///< Register unknown Pythia 8 particle to CepGen

  private:
    bool fillEventBlock();  ///< Copy all active central system particles into the Pythia 8 event record

    Event& cepgen_event_;      ///< CepGen event content (NOT owning)
    Pythia8::Pythia& pythia_;  ///< Pythia 8 engine (NOT owning)

    std::unordered_map<int, int> cepgen_ids_;  ///< Pythia 8 event record position -> CepGen identifier
  };
}  // namespace cepgen::pythia8

#endif
