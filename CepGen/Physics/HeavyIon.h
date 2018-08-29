#ifndef CepGen_Physics_HeavyIon_h
#define CepGen_Physics_HeavyIon_h

#include <ostream>

namespace CepGen
{
  enum class PDG;
  struct HeavyIon {
    static inline HeavyIon Proton() { return HeavyIon{ 1, 1 }; }
    static inline HeavyIon Pb208() { return HeavyIon{ 208, 82 }; }
    static HeavyIon fromPDG( const PDG& pdg );
    operator PDG() const;
    inline operator bool() const { return Z > 1; }
    friend std::ostream& operator<<( std::ostream& os, const HeavyIon& hi );
    /// Mass number
    unsigned short A;
    /// Atomic number
    unsigned short Z;
  };
}

#endif
