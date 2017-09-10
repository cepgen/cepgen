#ifndef CepGen_Physics_StructureFunctions_h
#define CepGen_Physics_StructureFunctions_h

#include "Particle.h"

extern "C"
{
  extern void grv95lo_( float&, float&, float&, float&, float&, float&, float&, float& );
}

namespace CepGen
{
  /// Proton structure function to be used in the outgoing state description
  /// \note Values correspond to the LPAIR legacy steering card values
  enum StructureFunctionsType {
    Electron = 1,
    ElasticProton = 2,
    SuriYennie = 11,
    SuriYennieLowQ2 = 12,
    SzczurekUleshchenko = 15,
    FioreVal = 101,
    FioreSea = 102,
    Fiore = 103
  };
  /// Human-readable format of a structure function type
  std::ostream& operator<<( std::ostream& os, const StructureFunctionsType& sf );

  class StructureFunctions
  {
    public:
      StructureFunctions( double f1=0.0, double f2=0.0 ) : F1( f1 ), F2( f2 ) {}
      /// Fiore-Brasse proton structure functions (F.W Brasse et al., DESY 76/11 (1976),
      /// http://dx.doi.org/10.1016/0550-3213(76)90231-5)
      /// \param[in] q2 Squared 4-momentum transfer
      /// \param[in] xbj Bjorken's x
      /// \cite Brasse1976413
      static StructureFunctions FioreBrasse( double q2, double xbj );
      static StructureFunctions SzczurekUleshchenko( double q2, double xbj );
      double F1, F2;
  };
  /// Human-readable format of a structure function object
  std::ostream& operator<<( std::ostream& os, const StructureFunctions& sf );
}

#endif
