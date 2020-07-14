#ifndef CepGen_Modules_EventModifierFactory_h
#define CepGen_Modules_EventModifierFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add a generic event modification module definition to the factory
#define REGISTER_MODIFIER( name, obj ) \
  namespace cepgen { namespace hadr { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { EventModifierFactory::get().registerModule<obj>( name ); } }; \
    static BUILDERNM( obj ) gEveMod ## obj; \
  } }

namespace cepgen
{
  class EventModifier;
  /// A event modifier algorithms factory
  typedef ModuleFactory<EventModifier> EventModifierFactory;
}

#endif

