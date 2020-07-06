#ifndef CepGen_Modules_FunctionalFactory_h
#define CepGen_Modules_FunctionalFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add a generic functional object builder definition
#define REGISTER_FUNCTIONAL( name, obj ) \
  namespace cepgen { namespace utils { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { FunctionalFactory::get().registerModule<obj>( name ); } }; \
    static BUILDERNM( obj ) gFunct ## obj; \
  } }

namespace cepgen
{
  namespace utils
  {
    class Functional;
    /// A functional objects factory
    typedef ModuleFactory<Functional> FunctionalFactory;
  }
}

#endif

