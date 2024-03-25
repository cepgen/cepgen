/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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

#include "CepGen/Modules/DocumentationGeneratorFactory.h"
#include "CepGen/Utils/DocumentationGenerator.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    /// Text documentation generator object
    /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
    /// \date Nov 2021
    class TextDocumentationGenerator final : public DocumentationGenerator {
    public:
      explicit TextDocumentationGenerator(const ParametersList& params)
          : DocumentationGenerator(params), dump_params_(steer<bool>("dumpParameters")) {}

      static ParametersDescription description() {
        auto desc = DocumentationGenerator::description();
        desc.setDescription("Bare text documentation generator");
        desc.add<bool>("dumpParameters", false)
            .setDescription("dump the parameters list along with their parameters description?");
        return desc;
      }

      std::string describe() override {
        std::ostringstream os;
        const auto separator = std::string(80, '=');
        for (const auto& cat : categories_) {
          if (cat.second.modules.empty())
            continue;
          os << "\n"
             << colourise(
                    separator + "\n" + cat.second.title + " modules" + "\n" + separator, Colour::green, Modifier::bold)
             << "\n";
          for (const auto& mod : cat.second.modules) {
            os << "\n"
               << colourise(toString(mod.first), Colour::none, Modifier::underline | Modifier::bold) << " module:\n\n"
               << mod.second.describe();
            if (dump_params_)
              os << "\n\tParametersList object:\n\t\t" << mod.second.parameters();
            os << "\n";
          }
        }
        return os.str();
      }

    private:
      const bool dump_params_;
    };
  }  // namespace utils
}  // namespace cepgen
using cepgen::utils::TextDocumentationGenerator;
REGISTER_DOCUMENTATION_GENERATOR("text", TextDocumentationGenerator);
