/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

/// Add a generic collinear parton flux evaluator builder definition
#define REGISTER_COLLINEAR_FLUX(name, obj)                                           \
  namespace cepgen {                                                                 \
    struct BUILDER_NAME(obj) {                                                       \
      BUILDER_NAME(obj)() { CollinearFluxFactory::get().registerModule<obj>(name); } \
    };                                                                               \
    static const BUILDER_NAME(obj) gCollinearFlux##obj;                              \
  }                                                                                  \
  static_assert(true, "")
/// Add a generic KT-factorised flux evaluator builder definition
#define REGISTER_KT_FLUX(name, id, obj)                                                          \
  namespace cepgen {                                                                             \
    struct BUILDER_NAME(obj) {                                                                   \
      BUILDER_NAME(obj)() { KTFluxFactory::get().addIndex(id, name).registerModule<obj>(name); } \
    };                                                                                           \
    static const BUILDER_NAME(obj) gKTFlux##obj;                                                 \
  }                                                                                              \
  static_assert(true, "")

namespace cepgen {
  class CollinearFlux;
  class KTFlux;
  /// A collinear parton fluxes objects factory
  DEFINE_FACTORY(CollinearFluxFactory, CollinearFlux, "Collinear parton flux estimators factory");
  /// A KT-factorised parton fluxes objects factory
  DEFINE_FACTORY(KTFluxFactory, KTFlux, "KT-factorised flux estimators factory");

  /// A generic parton fluxes objects factory
  struct PartonFluxFactory {
    static PartonFluxFactory& get() {
      static PartonFluxFactory instance;
      return instance;
    }

    static ParametersDescription describeParameters(const std::string& name,
                                                    const ParametersList& params = ParametersList());
    static bool elastic(const ParametersList&);           ///< Is the beam modelling elastic?
    static long long partonPdgId(const ParametersList&);  ///< Type of parton exchanged
  };
}  // namespace cepgen

#endif
