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
#include <iostream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Version.h"
#include "CepGenAddOns/CTMLWrapper/DocumentationGenerator.h"

namespace cepgen {
  namespace utils {
    DocumentationGenerator::DocumentationGenerator(const ParametersList& params)
        : SteeredObject(params),
          output_filename_(steer<std::string>("output")),
          container_(CTML::Node("div.container-fluid")) {
      doc_.AppendNodeToHead(CTML::Node("title", "CepGen v" + version::tag + " modules documentation"));
      if (steer<bool>("useBS")) {
        doc_.AppendNodeToHead(
            CTML::Node("link")
                .SetAttribute("rel", "stylesheet")
                .SetAttribute("href", "https://cdn.jsdelivr.net/npm/bootstrap@4.3.1/dist/css/bootstrap.min.css")
                .SetAttribute("integrity", "sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T")
                .SetAttribute("crossorigin", "anonymous"));
        doc_.AppendNodeToHead(CTML::Node("meta")
                                  .SetAttribute("name", "viewport")
                                  .SetAttribute("content", "width=device-width, initial-scale=1"));
      }
      container_.AppendChild(CTML::Node("h1", "Modules documentation"));
      container_.AppendChild(CTML::Node("div")
                                 .AppendText("CepGen version ")
                                 .AppendChild(CTML::Node("mark", version::tag))
                                 .AppendChild(CTML::Node("br").UseClosingTag(false))
                                 .AppendText("Git hash/branch: ")
                                 .AppendChild(CTML::Node("code", version::extended)));
    }

    DocumentationGenerator::~DocumentationGenerator() {
      doc_.AppendNodeToBody(container_);
      std::ostream* out{nullptr};
      if (!output_filename_.empty())
        out = new std::ofstream(output_filename_);
      else
        out = &std::cout;
      (*out) << doc_.ToString();
      if (!output_filename_.empty()) {
        delete out;
        CG_INFO("DocumentationGenerator") << "Documentation written in \"" << output_filename_ << "\".";
      }
    }

    CTML::Node DocumentationGenerator::moduleDescription(const ParametersDescription& desc) {
      CTML::Node out("div.module");
      if (desc.empty())
        return out;
      out.AppendChild(CTML::Node("b", desc.parameters().getString(ParametersList::MODULE_NAME)))
          .AppendText(" " + desc.description());
      try {
        CTML::Node items("ul");
        for (const auto& key : desc.parameters().keys(false)) {
          const auto& subdesc = desc.get(key);
          const auto type = subdesc.type();
          CTML::Node item("li.key");
          item.AppendChild(CTML::Node("u.key", key));
          if (type == ParametersDescription::Type::Value) {
            if (!subdesc.description().empty())
              item.AppendChild(CTML::Node("i", " " + subdesc.description()));
            if (!desc.parameters().getString(key).empty())
              item.AppendText(" ").AppendChild(
                  CTML::Node("span.text-muted")
                      .AppendText("(default value: ")
                      .AppendChild(CTML::Node("code", desc.parameters().getString(key, false)))
                      .AppendText(")"));
          } else
            item.AppendChild(moduleDescription(subdesc));
          items.AppendChild(item);
        }
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
      desc.add<std::string>("output", "index.html").setDescription("output path for the generated HTML file");
      desc.add<bool>("useBS", true).setDescription("use the Bootstrap CDN to prettify this output?");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen
