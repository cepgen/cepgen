/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2025  Laurent Forthomme
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

using namespace cepgen;
using namespace cepgen::utils;
using namespace std::string_literals;

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
    desc.add("modulesOnly", false).setDescription("only list the module names (for a category)?");
    desc.add("camelCaseModulesNames", false).setDescription("write modules in camel case?");
    desc.add("light", false).setDescription("lightweight module description (without parameters)");
    desc.add("dumpParameters", false).setDescription("dump parameters list along with their parameters description?");
    return desc;
  }

  std::string describe() override {
    std::ostringstream os;
    const auto separator = std::string(80, '-');
    const auto light = steer<bool>("light"), camel_case = steer<bool>("camelCaseModulesNames");
    std::vector<std::string> modules_names;
    for (const auto& [name, category] : categories_) {
      if (category.modules.empty())
        continue;
      os << colourise("\n" + separator + "\n" + category.title, Colour::green, Modifier::bold);
      if (!light)
        os << "\n";
      for (const auto& [module_name, description] : category.modules) {
        modules_names.emplace_back(camel_case ? toCamelCase(module_name) : module_name);
        if (light) {
          os << "\n"
             << (category.modules_indices.count(module_name) > 0
                     ? "#"s + std::to_string(category.modules_indices.at(module_name)) + ": "
                     : ""s)
             << colourise(module_name, Colour::cyan, Modifier::underline | Modifier::bold) << ": "
             << description.description() << (description.empty() ? " (*)" : "");
        } else {
          os << "\n";
          if (category.modules_indices.count(module_name) > 0)
            os << "#" << category.modules_indices.at(module_name) << ": ";
          os << description.describe();
          if (dump_params_)
            os << "\n\tParametersList object:\n\t\t" << description.parameters();
          os << "\n";
        }
      }
    }
    if (steer<bool>("modulesOnly"))
      return repr(modules_names, ";");
    return os.str();
  }

private:
  const bool dump_params_;
};
REGISTER_DOCUMENTATION_GENERATOR("text", TextDocumentationGenerator);
