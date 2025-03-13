/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#include <CTML/ctml.hpp>
#include <fstream>
#include <iostream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DocumentationGeneratorFactory.h"
#include "CepGen/Utils/DocumentationGenerator.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

namespace cepgen::ctml {
  /// CTML documentation generator object
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Apr 2022
  class DocumentationGenerator final : public utils::DocumentationGenerator {
  public:
    explicit DocumentationGenerator(const ParametersList& params)
        : utils::DocumentationGenerator(params),
          bare_(steer<bool>("bare")),
          container_(CTML::Node("div.container-fluid")) {}

    static ParametersDescription description() {
      auto desc = cepgen::utils::DocumentationGenerator::description();
      desc.setDescription("CTML HTML document generator helper");
      desc.add<std::string>("output", "index.html").setDescription("output path for the generated HTML file");
      desc.add<std::string>("pageTitle", "Modules documentation")
          .setDescription("documentation page upper level title");
      desc.add<bool>("useBS", true).setDescription("use the Bootstrap CDN to prettify this output?");
      desc.add<bool>("showGit", false).setDescription("print out the git hash/branch in the output?");
      desc.add<bool>("bare", false).setDescription("generate a bare version (without <html>/<head>/<body> attributes)");
      return desc;
    }

    std::string describe() override {
      doc_.AppendNodeToHead(CTML::Node("title", "CepGen v" + version::tag + " modules documentation"));
      if (!bare_ && steer<bool>("useBS")) {
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
      if (const auto page_title = steer<std::string>("pageTitle"); !page_title.empty())
        container_.AppendChild(CTML::Node("h1", page_title));
      auto header = CTML::Node("div").AppendText("CepGen version ").AppendChild(CTML::Node("mark", version::tag));
      if (steer<bool>("showGit"))
        header.AppendChild(CTML::Node("br").UseClosingTag(false))
            .AppendText("Git hash/branch: ")
            .AppendChild(CTML::Node("code", version::extended));
      header.AppendChild(CTML::Node("br").UseClosingTag(false))
          .AppendText("Documentation last generated on " + utils::timeAs("%B %d, %Y"));
      container_.AppendChild(header);
      for (const auto& cat : categories_) {
        container_.AppendChild(CTML::Node("a")
                                   .SetAttribute("name", cat.first)
                                   .AppendChild(CTML::Node("h2", cat.second.title).SetAttribute("id", cat.first)));
        CTML::Node mods("p");
        for (const auto& mod : cat.second.modules)
          mods.AppendChild(CTML::Node("a")
                               .SetAttribute("name", utils::toString(cat.first) + mod.first)
                               .AppendChild(CTML::Node("span").AppendChild(
                                   moduleDescription(mod.second,
                                                     cat.second.modules_indices.count(mod.first) > 0
                                                         ? cat.second.modules_indices.at(mod.first)
                                                         : -1)
                                       .SetAttribute("id", cat.first + mod.first))));
        container_.AppendChild(mods);
      }
      doc_.AppendNodeToBody(container_);
      std::ostringstream out;
      if (bare_)
        out << container_.ToString();
      else
        out << doc_.ToString();
      return out.str();
    }

  private:
    inline static CTML::Node moduleDescription(const ParametersDescription& desc, int index = -1) {
      CTML::Node out("div.module");
      if (desc.empty())
        return out;
      CTML::Node node_summary("summary");
      node_summary.AppendChild(CTML::Node("b", desc.parameters().getNameString()));
      if (index > 0)
        node_summary.AppendText(" (index ").AppendChild(CTML::Node("code", std::to_string(index))).AppendText(")");
      CTML::Node mod_summary(node_summary);
      CTML::Node mod_details("details");
      CTML::Node mod_params_list(CTML::Node("p").AppendText("List of parameters:"));
      const auto desc_type = desc.type();
      if (desc_type == ParametersDescription::Type::ParametersVector)
        mod_summary.AppendChild(CTML::Node("b", "Children attributes"));
      else if (desc_type == ParametersDescription::Type::Parameters)
        mod_summary.AppendChild(CTML::Node("b", "parameters list"));
      else if (desc_type == ParametersDescription::Type::Value) {
      } else {
        mod_summary.AppendText(" " + desc.description());
      }
      mod_details.AppendChild(mod_summary);
      try {
        CTML::Node items("ul");
        for (const auto& key : desc.parameters().keys(false)) {
          const auto& subdesc = desc.get(key);
          const auto subdesc_type = subdesc.type();
          CTML::Node item("li.key");
          item.AppendChild(CTML::Node("u.key", key));
          if (subdesc_type == ParametersDescription::Type::Value) {
            if (!subdesc.description().empty())
              item.AppendChild(CTML::Node("i", " " + subdesc.description()));
            if (!desc.parameters().getString(key).empty())
              item.AppendText(" ").AppendChild(
                  CTML::Node("span.text-muted")
                      .AppendText("(default value: ")
                      .AppendChild(CTML::Node("code", desc.parameters().getString(key, false)))
                      .AppendText(")"));
            if (const auto& allowed_vals = desc.get(key).allowedValues(); !allowed_vals.empty()) {
              item.AppendText(". Allowed values:");
              CTML::Node itparams("ul");
              for (const auto& it : allowed_vals.allowed()) {
                CTML::Node val("li");
                val.AppendChild(CTML::Node("code", it.first));
                if (!it.second.empty())
                  val.AppendText(" (" + it.second + ")");
                itparams.AppendChild(val);
              }
              item.AppendChild(itparams);
            }
          } else if (subdesc_type == ParametersDescription::Type::ParametersVector) {
            item.AppendText(" vector of parameters");
            if (!subdesc.description().empty())
              item.AppendText(" defining a ").AppendChild(CTML::Node("i", subdesc.description()));
            item.AppendChild(moduleDescription(subdesc));
            if (const auto& vparams = desc.parameters().get<std::vector<ParametersList> >(key); !vparams.empty()) {
              CTML::Node itparams("ol");
              for (const auto& it : vparams)
                itparams.AppendChild(CTML::Node("li").AppendChild(moduleDescription(ParametersDescription(it))));
              item.AppendChild(
                  CTML::Node("details")
                      .AppendChild(CTML::Node("summary").AppendChild(CTML::Node("b", "Default vector content")))
                      .AppendChild(CTML::Node("p").AppendChild(itparams.SetAttribute("start", "0"))));
            }
          } else
            item.AppendChild(CTML::Node("i", " " + subdesc.description())).AppendChild(moduleDescription(subdesc));
          items.AppendChild(item);
        }
        if (!items.GetChildren().empty()) {
          if (desc_type == ParametersDescription::Type::ParametersVector ||
              desc_type == ParametersDescription::Type::Parameters)
            mod_details.AppendChild(items);
          else
            mod_details.AppendChild(mod_params_list.AppendChild(items));
        }
        out.AppendChild(mod_details);
      } catch (const Exception& exc) {
        exc.dump();
      }
      return out;
    }

    const bool bare_;
    CTML::Document doc_;
    CTML::Node container_;
  };
}  // namespace cepgen::ctml
using DocumentationGeneratorCTML = cepgen::ctml::DocumentationGenerator;
REGISTER_DOCUMENTATION_GENERATOR("ctml", DocumentationGeneratorCTML);
