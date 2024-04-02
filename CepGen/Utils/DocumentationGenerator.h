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

#ifndef CepGen_Utils_DocumentationGenerator_h
#define CepGen_Utils_DocumentationGenerator_h

#include "CepGen/Modules/ModuleFactory.h"

namespace cepgen {
  namespace utils {
    /// Documentation generator object
    /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
    /// \date Mar 2024
    class DocumentationGenerator : public NamedModule<DocumentationGenerator> {
    public:
      explicit DocumentationGenerator(const ParametersList&);
      virtual ~DocumentationGenerator() = default;

      static ParametersDescription description();

      virtual std::string describe() = 0;

    protected:
      struct category_t {
        std::string name, title, description;
        std::map<std::string, ParametersDescription> modules{};
        std::map<std::string, int> modules_indices{};
      };
      std::vector<std::pair<std::string, category_t> > categories_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
