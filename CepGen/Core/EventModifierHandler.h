#ifndef CepGen_Core_EventModifierHandler_h
#define CepGen_Core_EventModifierHandler_h

#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/Core/EventModifier.h"

#define REGISTER_MODIFIER( name, obj ) \
  namespace cepgen { namespace hadr { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { EventModifierHandler::get().registerModule<obj>( name ); } }; \
    static BUILDERNM( obj ) g ## obj; \
  } }

namespace cepgen
{
  /// A event modifier algorithms factory
  typedef ModuleFactory<EventModifier> EventModifierHandler;
}

#endif

