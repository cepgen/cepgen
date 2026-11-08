/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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
#include "CepGenPython/ConfigWriter.h"

using namespace std::string_literals;

namespace cepgen::python {
  /// Python modules documentation generator
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Mar 2024
  class DocumentationGenerator final : public utils::DocumentationGenerator {
  public:
    explicit DocumentationGenerator(const ParametersList& params)
        : utils::DocumentationGenerator(params), writer_(params_) {}

    static ParametersDescription description() {
      auto desc = utils::DocumentationGenerator::description();
      desc.setDescription("Python modules documentation generator");
      desc += ConfigWriter::description();
      desc.add("filename", ""s).setDescription("Python output filename");
      desc.add("camelCaseModuleNames", true);
      return desc;
    }

    std::string describe() override {
      for (const auto& [name, category] : categories_)
        for (const auto& [module_name, description] : category.modules)
          writer_ << description;
      return writer_();
    }

  private:
    ConfigWriter writer_;
  };
}  // namespace cepgen::python
using PythonDocumentationGenerator = cepgen::python::DocumentationGenerator;
REGISTER_DOCUMENTATION_GENERATOR("python", PythonDocumentationGenerator);
