#ifndef CepGen_Physics_HeavyIon_h
#define CepGen_Physics_HeavyIon_h

#include <ostream>

namespace cepgen
{
  enum class PDG;
  /// Enumeration of chemical elements
  enum class Element
  {
    invalid = 0,
    H = 1, C = 6, O = 8,
    Al = 13, Cu = 29,
    Xe = 54, Au = 79, Pb = 82,
    U = 92
  };
  std::ostream& operator<<( std::ostream& os, const Element& elem );

  /// Heavy ion container (Z+A)
  struct HeavyIon
  {
    /// General constructor from mass and atomic number
    HeavyIon( unsigned short a, const Element& z ) : A( a ), Z( z ) {}
    /// Build from a custom PDG id
    HeavyIon( const PDG& pdg );
    /// Simple proton
    static inline HeavyIon proton() { return HeavyIon( 1, Element::H ); }
    /// Convert the HI into a custom PDG id
    operator PDG() const;
    /// Check the validity of the heavy ion
    inline operator bool() const {
      return Z > Element::invalid && A > 1; // skip the proton
    }
    /// Human-readable expression of the ion
    friend std::ostream& operator<<( std::ostream& os, const HeavyIon& hi );
    /// Mass number
    unsigned short A;
    /// Atomic number
    Element Z;
  };
}

#endif
