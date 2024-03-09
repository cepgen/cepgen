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
  enum class Element {
    invalid = -1,
    neutron = 0,
    H = 1,
    C = 6,
    O = 8,
    Al = 13,
    Cu = 29,
    Xe = 54,
    Au = 79,
    Pb = 82,
    U = 92
  };
  std::ostream& operator<<(std::ostream& os, const Element& elem);

  /// Heavy ion container (Z+A)
  struct HeavyIon {
    /// General constructor from mass and atomic number
    explicit HeavyIon(unsigned short, const Element&);

    bool operator==(const HeavyIon& oth) const { return Z == oth.Z && A == oth.A; }
    bool operator!=(const HeavyIon& oth) const { return !(*this == oth); }

    /// Neutrons mass, in GeV/c2
    double massN() const;
    /// Protons mass, in GeV/c2
    double massP() const;
    /// Total heavy ion mass, in GeV/c2
    double mass() const { return massN() + massP(); }
    /// Heavy ion radius, in m
    double radius() const;

    /// Build from a custom PDG id
    static HeavyIon fromPdgId(pdgid_t);
    /// Check if the PDG id is compatible with a HI
    static bool isHI(const spdgid_t&);
    /// Check if the particle properties are compatible with a HI
    static bool isHI(const ParticleProperties&);
    /// Mass of a heavy ion, in GeV/c\f$^2\f$
    /// \param hi Heavy ion type
    static double mass(const HeavyIon& hi) { return hi.mass(); }

    /// Simple proton
    static inline HeavyIon proton() { return HeavyIon(1, Element::H); }
    /// Simple neutron
    static inline HeavyIon neutron() { return HeavyIon(1, Element::neutron); }
    /// Standard gold
    static inline HeavyIon Au() { return HeavyIon(197, Element::Au); }
    /// Standard lead
    static inline HeavyIon Pb() { return HeavyIon(207, Element::Pb); }

    /// Convert the HI into a custom PDG id
    operator pdgid_t() const;
    /// Human-readable expression of the ion
    friend std::ostream& operator<<(std::ostream& os, const HeavyIon& hi);

    /// Atomic number
    Element Z{Element::invalid};
    /// Mass number
    unsigned short A{0};
  };
}  // namespace cepgen

#endif
