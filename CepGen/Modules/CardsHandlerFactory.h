/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

/// Add a cards handler definition to the list of handled parsers
#define REGISTER_CARD_HANDLER(name, obj)                                           \
  namespace cepgen {                                                               \
    namespace card {                                                               \
      struct BUILDERNM(obj) {                                                      \
        BUILDERNM(obj)() { CardsHandlerFactory::get().registerModule<obj>(name); } \
      };                                                                           \
      static const BUILDERNM(obj) gCard##obj;                                      \
    }                                                                              \
  }                                                                                \
  static_assert(true, "")

namespace cepgen {
  class RunParameters;
  namespace card {
    class Handler;
  }
  /// A cards handler base factory
  DEFINE_FACTORY(BaseCardsHandlerFactory, card::Handler, "Cards handlers factory");
  /// A cards handler factory
  struct CardsHandlerFactory : public BaseCardsHandlerFactory {
    using BaseCardsHandlerFactory::BaseCardsHandlerFactory;
    static CardsHandlerFactory& get();
    /// Build one instance of a cards handler
    /// \param[in] filename File path to retrieve
    std::unique_ptr<card::Handler> buildFromFilename(const std::string& filename) const;
  };
}  // namespace cepgen

#endif
