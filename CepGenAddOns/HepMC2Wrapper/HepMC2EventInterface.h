/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#ifndef CepGenAddOns_HepMC2Wrapper_HepMC2EventInterface_h
#define CepGenAddOns_HepMC2Wrapper_HepMC2EventInterface_h

#include <HepMC/GenEvent.h>

#include <memory>
#include <unordered_map>

namespace cepgen {
  class Event;
}

namespace HepMC {
  /// Interfacing between CepGen and HepMC2 event definitions
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jul 2019
  class CepGenEvent : public GenEvent {
  public:
    /// Construct an event interface from a CepGen Event object
    CepGenEvent(const cepgen::Event& ev);
    /// Extract a CepGen Event object from a HepMC2 GenEvent object
    operator cepgen::Event() const;

  private:
    std::unordered_map<unsigned short, GenParticle*> assoc_map_;
  };
}  // namespace HepMC
#endif
