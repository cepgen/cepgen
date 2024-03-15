/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#ifndef CepGenAddOns_CTMLWrapper_DocumentationGenerator_h
#define CepGenAddOns_CTMLWrapper_DocumentationGenerator_h

#include <CTML/ctml.hpp>
#include <sstream>

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Modules/ModuleFactory.h"

namespace cepgen {
  namespace utils {
    /**
     * \brief CTML documentation generator object
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Apr 2022
     */
    class DocumentationGenerator final : public SteeredObject<DocumentationGenerator> {
    public:
      explicit DocumentationGenerator(const ParametersList&);
      ~DocumentationGenerator();

      static ParametersDescription description();

      template <typename T, typename I>
      inline DocumentationGenerator& document(const std::string& type,
                                              const std::string& title,
                                              const ModuleFactory<T, I>& factory) {
        container_.AppendChild(CTML::Node("a").SetAttribute("name", type)).AppendChild(CTML::Node("h2", title));
        CTML::Node mods("p");
        for (const auto& mod : factory.modules()) {
          std::ostringstream os;
          os << type << "-" << mod;
          mods.AppendChild(CTML::Node("a").SetAttribute("name", os.str()))
              .AppendChild(CTML::Node("span").AppendChild(moduleDescription(factory.describeParameters(mod))));
        }
        container_.AppendChild(mods);
        return *this;
      }

    private:
      static CTML::Node moduleDescription(const ParametersDescription&);

      const std::string output_filename_;
      const bool bare_, show_git_;
      CTML::Document doc_;
      CTML::Node container_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
