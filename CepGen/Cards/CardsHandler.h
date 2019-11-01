#ifndef CepGen_Cards_CardsHandler_h
#define CepGen_Cards_CardsHandler_h

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/ModuleFactory.h"

/** \file */

/// Add a cards handler definition to the list of handled parsers
#define REGISTER_CARD_HANDLER( name, obj ) \
  namespace cepgen { namespace card { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { CardsHandler::get().registerModule<obj>( name ); } }; \
    static BUILDERNM( obj ) g ## obj; \
  } }

namespace cepgen
{
  namespace card
  {
    /// A cards handler factory
    typedef ModuleFactory<Handler> CardsHandler;
  }
}

#endif

