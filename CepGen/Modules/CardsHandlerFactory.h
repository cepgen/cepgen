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

#ifndef CepGen_Modules_CardsHandlerFactory_h
#define CepGen_Modules_CardsHandlerFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add a cards handler definition to the list of handled parsers
#define REGISTER_CARD_HANDLER(name, obj)                                           \
  namespace cepgen {                                                               \
    namespace card {                                                               \
      struct BUILDERNM(obj) {                                                      \
        BUILDERNM(obj)() { CardsHandlerFactory::get().registerModule<obj>(name); } \
      };                                                                           \
      static const BUILDERNM(obj) gCard##obj;                                      \
    }                                                                              \
  }

namespace cepgen {
  namespace card {
    class Handler;
    /// A cards handler factory
    DEFINE_FACTORY_STR(CardsHandlerFactory, Handler, "Cards handlers factory");
  }  // namespace card
}  // namespace cepgen

#endif
