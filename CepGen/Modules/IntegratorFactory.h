#ifndef CepGen_Modules_IntegratorFactory_h
#define CepGen_Modules_IntegratorFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add a generic process definition to the list of handled processes
#define REGISTER_INTEGRATOR(name, obj)                                         \
  namespace cepgen {                                                           \
    struct BUILDERNM(obj) {                                                    \
      BUILDERNM(obj)() { IntegratorFactory::get().registerModule<obj>(name); } \
    };                                                                         \
    static const BUILDERNM(obj) gIntegr##obj;                                  \
  }

namespace cepgen {
  class Integrator;
  /// An integration algorithms factory
  DEFINE_FACTORY_STR(IntegratorFactory, Integrator, "Integrator factory");
}  // namespace cepgen

#endif
