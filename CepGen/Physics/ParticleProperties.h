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

#ifndef CepGen_Physics_ParticleProperties_h
#define CepGen_Physics_ParticleProperties_h

#include <iosfwd>
#include <string>

namespace cepgen {
  /// Alias for the integer-like particle PDG id
  typedef unsigned long pdgid_t;
  /// A collection of physics constants associated to a single particle
  struct ParticleProperties {
    pdgid_t pdgid;            ///< PDG identifier
    std::string name;         ///< Particle name
    std::string description;  ///< Human-readable name
    short colours;            ///< Colour factor
    double mass;              ///< Mass, in GeV/c\f$^2\f$
    double width;             ///< Decay width, in GeV/c\f$^2\f$
    short charge;             ///< Electric charge, in \f$e\f$/3
    bool fermion;             ///< Is the particle a fermion?
    friend std::ostream& operator<<(std::ostream&, const ParticleProperties&);
  };
}  // namespace cepgen

#endif
