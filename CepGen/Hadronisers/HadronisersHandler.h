#ifndef CepGen_Hadronisers_HadronisersHandler_h
#define CepGen_Hadronisers_HadronisersHandler_h

#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#define BUILDERNM( obj ) obj ## Builder
#define STRINGIFY( name ) #name
#define REGISTER_HADRONISER( name, obj ) \
  namespace cepgen { namespace hadr { \
    struct BUILDERNM( name ) { \
      BUILDERNM( name )() { HadronisersHandler::get().registerModule<obj>( STRINGIFY( name ) ); } }; \
    static BUILDERNM( name ) g ## name; \
  } }

namespace cepgen
{
  namespace hadr
  {
    /// A hadroniser modules factory
    typedef ModuleFactory<GenericHadroniser> HadronisersHandler;
  }
}

#endif

