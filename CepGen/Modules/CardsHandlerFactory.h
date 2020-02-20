#ifndef CepGen_Modules_CardsHandlerFactory_h
#define CepGen_Modules_CardsHandlerFactory_h

#include "CepGen/Core/ModuleFactory.h"

/** \file */

/// Add a cards handler definition to the list of handled parsers
#define REGISTER_CARD_HANDLER( name, obj ) \
  namespace cepgen { namespace card { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { CardsHandlerFactory::get().registerModule<obj>( name ); } }; \
    static BUILDERNM( obj ) gCard ## obj; \
  } }

namespace cepgen
{
  namespace card
  {
    class Handler;
    /// A cards handler factory
    typedef ModuleFactory<Handler> CardsHandlerFactory;
  }
}

#endif

