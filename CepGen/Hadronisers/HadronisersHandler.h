#ifndef CepGen_Hadronisers_HadronisersHandler_h
#define CepGen_Hadronisers_HadronisersHandler_h

#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#define REGISTER_HADRONISER( name, obj ) \
  namespace cepgen { namespace hadr { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { HadronisersHandler::get().registerModule<obj>( name ); } }; \
    static BUILDERNM( obj ) g ## obj; \
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

