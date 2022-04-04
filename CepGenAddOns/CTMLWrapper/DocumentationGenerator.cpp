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

#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Version.h"
#include "CepGenAddOns/CTMLWrapper/DocumentationGenerator.h"

namespace cepgen {
  namespace utils {
    DocumentationGenerator::DocumentationGenerator(const ParametersList& params)
        : SteeredObject(params), output_filename_(steer<std::string>("output")) {
      doc_.AppendNodeToHead(CTML::Node("title", "CepGen v" + version::tag + " modules documentation"));
      doc_.AppendNodeToBody(CTML::Node("h1", "Modules documentation"));
      doc_.AppendNodeToBody(CTML::Node("div", "CepGen version " + version::tag));
    }

    DocumentationGenerator::~DocumentationGenerator() {
      if (!output_filename_.empty()) {
        std::ofstream out(output_filename_);
        out << doc_.ToString();
        out.close();
        CG_INFO("DocumentationGenerator") << "Documentation written in \"" << output_filename_ << "\".";
      }
    }

    CTML::Node DocumentationGenerator::moduleDescription(const ParametersDescription& desc) {
      CTML::Node out("div");
      if (desc.empty())
        return out;
      out.AppendChild(CTML::Node("b", desc.parameters().getString(ParametersList::MODULE_NAME)))
          .AppendText(" " + desc.description());
      try {
        CTML::Node items("ul");
        for (const auto& key : desc.parameters().keys(false))
          items.AppendChild(CTML::Node("li", key).AppendChild(moduleDescription(desc.get(key))));
        if (!items.GetChildren().empty())
          out.AppendChild(items);
      } catch (const Exception& exc) {
        exc.dump();
      }
      return out;
    }

    ParametersDescription DocumentationGenerator::description() {
      auto desc = ParametersDescription();
      desc.setDescription("CTML HTML document generator helper");
      desc.add<std::string>("output", "index.html");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen
