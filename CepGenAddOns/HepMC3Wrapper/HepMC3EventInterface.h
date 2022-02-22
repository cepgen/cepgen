/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#ifndef CepGenAddOns_HepMC3Wrapper_HepMC3EventInterface_h
#define CepGenAddOns_HepMC3Wrapper_HepMC3EventInterface_h

#include <HepMC3/GenEvent.h>

#include <memory>
#include <unordered_map>

namespace cepgen {
  class Event;
}

namespace HepMC3 {
  /// Interfacing between CepGen and HepMC event definitions
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jul 2019
  class CepGenEvent : public GenEvent {
  public:
    /// Construct an event interface from a CepGen Event object
    CepGenEvent(const cepgen::Event& ev);

  private:
    std::unordered_map<unsigned short, std::shared_ptr<GenParticle> > assoc_map_;
  };
}  // namespace HepMC3
#endif
