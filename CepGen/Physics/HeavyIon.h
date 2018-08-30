#ifndef CepGen_Physics_HeavyIon_h
#define CepGen_Physics_HeavyIon_h

#include <ostream>

namespace CepGen
{
  enum class PDG;
  enum class Element
  {
    H = 1, C = 6, O = 8,
    Al = 13, Cu = 29,
    Xe = 54, Au = 79, Pb = 82,
    U = 92
  };
  std::ostream& operator<<( std::ostream& os, const Element& elem );

  /// Heavy ion container (Z+A)
  struct HeavyIon
  {
    static inline HeavyIon proton() { return HeavyIon{ 1, Element::H }; }
    operator PDG() const;
    inline operator bool() const { return Z > Element::H; }
    friend std::ostream& operator<<( std::ostream& os, const HeavyIon& hi );
    /// Mass number
    unsigned short A;
    /// Atomic number
    Element Z;
  };
}

#endif
