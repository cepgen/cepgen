/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#ifndef CepGen_Modules_DocumentationGeneratorFactory_h
#define CepGen_Modules_DocumentationGeneratorFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a documentation generator to the list of handled modules
#define REGISTER_DOCUMENTATION_GENERATOR(name, obj)                                        \
  namespace cepgen {                                                                       \
    struct BUILDERNM(obj) {                                                                \
      BUILDERNM(obj)() { DocumentationGeneratorFactory::get().registerModule<obj>(name); } \
    };                                                                                     \
    static const BUILDERNM(obj) gDogGen##obj;                                              \
  }                                                                                        \
  static_assert(true, "")

namespace cepgen {
  namespace utils {
    class DocumentationGenerator;
  }
  /// A documentation generator factory
  DEFINE_FACTORY(DocumentationGeneratorFactory, utils::DocumentationGenerator, "Documentation generator factory");
}  // namespace cepgen

#endif
