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

#include "CepGen/Modules/DocumentationGeneratorFactory.h"
#include "CepGen/Utils/DocumentationGenerator.h"
#include "CepGenAddOns/PythonWrapper/ConfigWriter.h"

using namespace cepgen;

/// Python modules documentation generator
/// \author Laurent Forthomme <laurent.forthomme@cern.ch>
/// \date Mar 2024
class PythonDocumentationGenerator final : public cepgen::utils::DocumentationGenerator {
public:
  explicit PythonDocumentationGenerator(const ParametersList& params) : DocumentationGenerator(params) {}

  static ParametersDescription description() {
    auto desc = DocumentationGenerator::description();
    desc.setDescription("Python modules documentation generator");
    desc.add<std::string>("filename", "output.py").setDescription("Python output filename");
    return desc;
  }

  std::string describe() override {
    python::ConfigWriter writer(params_);
    for (const auto& cat : categories_) {
      if (cat.second.modules.empty())
        continue;
      for (const auto& mod : cat.second.modules)
        writer << mod.second;
    }
    return "";
  }
};
REGISTER_DOCUMENTATION_GENERATOR("python", PythonDocumentationGenerator);
