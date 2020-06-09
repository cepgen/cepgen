#ifndef CepGen_Physics_AlphaS_h
#define CepGen_Physics_AlphaS_h

#include "CepGen/Core/ModuleFactory.h"

#define REGISTER_ALPHAS_MODULE( name, obj ) \
  namespace cepgen { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { AlphaSFactory::get().registerModule<obj>( name ); } }; \
    static BUILDERNM( obj ) g ## obj; \
  }

namespace cepgen
{
  /// A generic \f$\alpha_S(Q^2)\f$ evaluation algoritm
  class AlphaS
  {
    public:
      AlphaS() = default;
      /// Compute \f$\alpha_S\f$ for a given \f$Q^2\f$
      virtual double operator()( double q ) const = 0;
  };
  /// An alpha(S) evolution algorithms factory
  typedef ModuleFactory<AlphaS> AlphaSFactory;
}

#endif
