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

#ifndef CepGen_Physics_HeavyIon_h
#define CepGen_Physics_HeavyIon_h

#include <ostream>

#include "CepGen/Physics/ParticleProperties.h"

namespace cepgen {
  /// Enumeration of chemical elements
  enum class Element { invalid = 0, H = 1, C = 6, O = 8, Al = 13, Cu = 29, Xe = 54, Au = 79, Pb = 82, U = 92 };
  std::ostream& operator<<(std::ostream& os, const Element& elem);

  /// Heavy ion container (Z+A)
  struct HeavyIon {
    /// General constructor from mass and atomic number
    HeavyIon(unsigned short a, const Element& z) : Z(z), A(a) {}
    /// Build from a custom PDG id
    HeavyIon(pdgid_t pdg);
    /// Check if the PDG id is compatible with a HI
    static bool isHI(const pdgid_t&);
    /// Simple proton
    static inline HeavyIon proton() { return HeavyIon(1, Element::H); }
    /// Mass of a heavy ion, in GeV/c\f$^2\f$
    /// \param hi Heavy ion type
    static double mass(const HeavyIon& hi);
    /// Convert the HI into a custom PDG id
    operator pdgid_t() const;
    /// Check the validity of the heavy ion
    operator bool() const;
    /// Human-readable expression of the ion
    friend std::ostream& operator<<(std::ostream& os, const HeavyIon& hi);
    /// Atomic number
    Element Z;
    /// Mass number
    unsigned short A;
  };
}  // namespace cepgen

#endif
