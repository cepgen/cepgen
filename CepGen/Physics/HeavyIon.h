/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2025  Laurent Forthomme
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

#include "CepGen/Physics/ParticleProperties.h"

namespace cepgen {
  /// Enumeration of chemical elements
  enum class Element {
    invalid = -1,
    neutron = 0,
    H = 1,    ///< hydrogen
    C = 6,    ///< carbon
    O = 8,    ///< oxygen
    Al = 13,  ///< aluminium
    Cu = 29,  ///< copper
    Xe = 54,  ///< xenon
    Au = 79,  ///< gold
    Pb = 82,  ///< lead
    U = 92    ///< uranium
  };
  std::ostream& operator<<(std::ostream& os, const Element& elem);

  /// Heavy ion container (Z+A)
  struct HeavyIon {
    explicit HeavyIon(unsigned short, const Element&);  ///< General constructor from mass and atomic number

    bool operator==(const HeavyIon& oth) const { return Z == oth.Z && A == oth.A; }
    bool operator!=(const HeavyIon& oth) const { return !(*this == oth); }

    double massN() const;                                     ///< Mass of all neutrons in HI, in GeV/c2
    double massP() const;                                     ///< Mass of all protons in HI, in GeV/c2
    inline double mass() const { return massN() + massP(); }  ///< Total heavy ion mass, in GeV/c2
    double radius() const;                                    ///< Heavy ion radius, in m

    static HeavyIon fromPdgId(pdgid_t);           ///< Build a HI from a custom PDG id
    static bool isHI(const spdgid_t&);            ///< Check if the PDG id is compatible with a HI
    static bool isHI(const ParticleProperties&);  ///< Check if the particle properties are compatible with a HI

    /// Mass of a heavy ion, in GeV/c\f$^2\f$
    /// \param heavy_ion Heavy ion type
    static double mass(const HeavyIon& heavy_ion) { return heavy_ion.mass(); }

    static inline HeavyIon proton() { return HeavyIon(1, Element::H); }         ///< Simple proton
    static inline HeavyIon neutron() { return HeavyIon(1, Element::neutron); }  ///< Simple neutron
    static inline HeavyIon Au() { return HeavyIon(197, Element::Au); }          ///< Standard gold
    static inline HeavyIon Pb() { return HeavyIon(207, Element::Pb); }          ///< Standard lead

    operator pdgid_t() const;                                         ///< Convert the HI into a custom PDG id
    friend std::ostream& operator<<(std::ostream&, const HeavyIon&);  ///< Human-readable expression of the ion

    Element Z{Element::invalid};  ///< Atomic number
    unsigned short A{0};          ///< Mass number
  };
}  // namespace cepgen

#endif
