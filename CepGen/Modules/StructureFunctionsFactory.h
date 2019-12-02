#ifndef CepGen_Modules_StructureFunctionsFactory_h
#define CepGen_Modules_StructureFunctionsFactory_h

#include "CepGen/Core/ModuleFactory.h"

/// Add a structure functions definition to the list of handled parameterisation
#define REGISTER_STRFUN( id, obj ) \
  namespace cepgen { \
    struct BUILDERNM( id ) { \
      BUILDERNM( id )() { strfun::StructureFunctionsFactory::get().registerModule<obj>( (int)strfun::Type::id ); } }; \
    static BUILDERNM( id ) g ## id; \
  }
/// Add a structure functions definition (with its associated default parameters) to the list of handled parameterisation
#define REGISTER_STRFUN_PARAMS( id, obj, params ) \
  namespace cepgen { \
    struct BUILDERNM( id ) { \
      BUILDERNM( id )() { strfun::StructureFunctionsFactory::get().registerModule<obj>( (int)strfun::Type::id, params ); } }; \
    static BUILDERNM( id ) g ## id; \
  }

namespace cepgen
{
  namespace strfun
  {
    class Parameterisation;
    /// A structure functions parameterisations factory
    typedef ModuleFactory<Parameterisation,int> StructureFunctionsFactory;
  }
}

#endif
