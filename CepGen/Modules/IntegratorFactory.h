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
  /// A processes factory
  struct IntegratorFactory : public ModuleFactory<Integrator> {
    static std::string description() { return "Integrators factory"; }
  };
}  // namespace cepgen

#endif
