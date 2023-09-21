/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Modules_PartonFluxFactory_h
#define CepGen_Modules_PartonFluxFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add a generic KT-factorised flux evaluator builder definition
#define REGISTER_KT_FLUX(name, obj)                                        \
  namespace cepgen {                                                       \
    struct BUILDERNM(obj) {                                                \
      BUILDERNM(obj)() { KTFluxFactory::get().registerModule<obj>(name); } \
    };                                                                     \
    static const BUILDERNM(obj) gKTFlux##obj;                              \
  }                                                                        \
  static_assert(true, "")

namespace cepgen {
  class KTFlux;
  /// A KT-factorised parton fluxes objects factory
  DEFINE_FACTORY_STR(KTFluxFactory, KTFlux, "KT-factorised flux estimators factory");

  struct PartonFluxFactory {
    static PartonFluxFactory& get() {
      static PartonFluxFactory instance;
      return instance;
    }

    ParametersDescription describeParameters(const std::string& name,
                                             const ParametersList& params = ParametersList()) const;
  };
}  // namespace cepgen

#endif
