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
    HeavyIon( unsigned short a, const Element& z ) : Z( z ), A( a ) {}
    /// Build from a custom PDG id
    HeavyIon( const PDG& pdg );
    /// Simple proton
    static inline HeavyIon proton() { return HeavyIon( 1, Element::H ); }
    /// Convert the HI into a custom PDG id
    operator PDG() const;
    /// Check the validity of the heavy ion
    operator bool() const;
    /// Human-readable expression of the ion
    friend std::ostream& operator<<( std::ostream& os, const HeavyIon& hi );
    /// Atomic number
    Element Z;
    /// Mass number
    unsigned short A;
  };
  namespace particleproperties
  {
    /// Mass of a heavy ion, in GeV/c\f$^2\f$
    /// \param hi Heavy ion type
    double mass( const HeavyIon& hi );
  }
}

#endif
