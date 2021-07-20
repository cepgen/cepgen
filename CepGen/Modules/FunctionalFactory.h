#ifndef CepGen_Modules_FunctionalFactory_h
#define CepGen_Modules_FunctionalFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add a generic functional object builder definition
#define REGISTER_FUNCTIONAL(name, obj)                                           \
  namespace cepgen {                                                             \
    namespace utils {                                                            \
      struct BUILDERNM(obj) {                                                    \
        BUILDERNM(obj)() { FunctionalFactory::get().registerModule<obj>(name); } \
      };                                                                         \
      static const BUILDERNM(obj) gFunct##obj;                                   \
    }                                                                            \
  }

namespace cepgen {
  namespace utils {
    class Functional;
    /// A functional objects factory
    DEFINE_FACTORY(FunctionalFactory, Functional, "Functionals factory");
  }  // namespace utils
}  // namespace cepgen

#endif
