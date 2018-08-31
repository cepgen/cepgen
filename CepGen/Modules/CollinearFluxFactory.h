/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#ifndef CepGen_Modules_CollinearFluxFactory_h
#define CepGen_Modules_CollinearFluxFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a collinear flux definition to the list of handled parameterisation
#define REGISTER_COLLFLUX(id, obj)                                                                              \
  namespace cepgen {                                                                                            \
    namespace collflux {                                                                                        \
      struct BUILDERNM(id) {                                                                                    \
        BUILDERNM(id)() { collflux::CollinearFluxFactory::get().registerModule<obj>((int)collflux::Type::id); } \
      };                                                                                                        \
      static const BUILDERNM(id) gCollFlux##id;                                                                 \
    }                                                                                                           \
  }

namespace cepgen {
  namespace collflux {
    class Parameterisation;
  }
  /// A collinear flux parameterisations factory
  DEFINE_FACTORY_INT(CollinearFluxFactory, collflux::Parameterisation, "Collinear flux parameterisations factory");
}  // namespace cepgen

#endif
