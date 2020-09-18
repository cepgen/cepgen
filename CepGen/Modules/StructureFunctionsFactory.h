#ifndef CepGen_Modules_StructureFunctionsFactory_h
#define CepGen_Modules_StructureFunctionsFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a structure functions definition to the list of handled parameterisation
#define REGISTER_STRFUN( id, obj ) \
  namespace cepgen { \
    struct BUILDERNM( id ) { \
      BUILDERNM( id )() { strfun::StructureFunctionsFactory::get().registerModule<obj>( (int)strfun::Type::id ); } }; \
    static BUILDERNM( id ) gStrFun ## id; \
  }
/// Add a structure functions definition (with its associated default parameters) to the list of handled parameterisation
#define REGISTER_STRFUN_PARAMS( id, obj, params ) \
  namespace cepgen { \
    struct BUILDERNM( id ) { \
      BUILDERNM( id )() { strfun::StructureFunctionsFactory::get().registerModule<obj>( (int)strfun::Type::id, params ); } }; \
    static BUILDERNM( id ) gStrFun ## id; \
  }

/// Add a sigma ratio definition to the list of handled parameterisation
#define SRBUILDERNM( id ) Ratio ## id ## Builder
#define REGISTER_SIGRAT( id, obj ) \
  namespace cepgen { \
    struct SRBUILDERNM( id ) { \
      SRBUILDERNM( id )() { sigrat::SigmaRatiosFactory::get().registerModule<obj>( (int)sigrat::Type::id ); } }; \
    static SRBUILDERNM( id ) gSigRat ## id; \
  }

/// Add a form factors definition to the list of handled parameterisation
#define REGISTER_FF_MODEL( id, obj ) \
  namespace cepgen { \
    struct BUILDERNM( id ) { \
      BUILDERNM( id )() { formfac::FormFactorsFactory::get().registerModule<obj>( (int)formfac::Model::id ); } }; \
    static BUILDERNM( id ) gFF ## id; \
  }

namespace cepgen
{
  namespace strfun
  {
    class Parameterisation;
    /// A structure functions parameterisations factory
    typedef ModuleFactory<Parameterisation,int> StructureFunctionsFactory;
  }
  namespace sigrat
  {
    class Parameterisation;
    /// A sigma ratio parameterisations factory
    typedef ModuleFactory<Parameterisation,int> SigmaRatiosFactory;
  }
  namespace formfac
  {
    class Parameterisation;
    /// A form factors parameterisations factory
    typedef ModuleFactory<Parameterisation,int> FormFactorsFactory;
  }
}

#endif
