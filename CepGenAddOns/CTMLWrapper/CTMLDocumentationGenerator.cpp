/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

using namespace cepgen;
/// CTML documentation generator object
/// \author Laurent Forthomme <laurent.forthomme@cern.ch>
/// \date Apr 2022
class CTMLDocumentationGenerator final : public cepgen::utils::DocumentationGenerator {
public:
  explicit CTMLDocumentationGenerator(const ParametersList& params)
      : DocumentationGenerator(params),
        page_title_(steer<std::string>("pageTitle")),
        bare_(steer<bool>("bare")),
        show_git_(steer<bool>("showGit")),
        container_(CTML::Node("div.container-fluid")) {
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
    if (!page_title_.empty())
      container_.AppendChild(CTML::Node("h1", page_title_));
    auto header = CTML::Node("div").AppendText("CepGen version ").AppendChild(CTML::Node("mark", version::tag));
    if (show_git_)
      header.AppendChild(CTML::Node("br").UseClosingTag(false))
          .AppendText("Git hash/branch: ")
          .AppendChild(CTML::Node("code", version::extended));
    header.AppendChild(CTML::Node("br").UseClosingTag(false))
        .AppendText("Documentation last generated on " + utils::timeAs("%B %d, %Y"));
    container_.AppendChild(header);
  }

  static ParametersDescription description() {
    auto desc = DocumentationGenerator::description();
    desc.setDescription("CTML HTML document generator helper");
    desc.add<std::string>("output", "index.html").setDescription("output path for the generated HTML file");
    desc.add<std::string>("pageTitle", "Modules documentation").setDescription("documentation page upper level title");
    desc.add<bool>("useBS", true).setDescription("use the Bootstrap CDN to prettify this output?");
    desc.add<bool>("showGit", false).setDescription("print out the git hash/branch in the output?");
    desc.add<bool>("bare", false).setDescription("generate a bare version (without <html>/<head>/<body> attributes)");
    return desc;
  }

  void initialise() override {
    DocumentationGenerator::initialise();
    for (const auto& cat : categories_) {
      container_.AppendChild(CTML::Node("a").SetAttribute("name", cat.first))
          .AppendChild(CTML::Node("h2", cat.second.title));
      CTML::Node mods("p");
      for (const auto& mod : cat.second.modules)
        mods.AppendChild(CTML::Node("a").SetAttribute("name", utils::toString(cat.first) + mod.first))
            .AppendChild(CTML::Node("span").AppendChild(moduleDescription(mod.second)));
      container_.AppendChild(mods);
    }
  }
  std::string describe() override {
    doc_.AppendNodeToBody(container_);
    std::ostringstream out;
    if (bare_)
      out << container_.ToString();
    else
      out << doc_.ToString();
    return out.str();
  }

private:
  inline static CTML::Node moduleDescription(const ParametersDescription& desc) {
    CTML::Node out("div.module");
    if (desc.empty())
      return out;
    CTML::Node mod_summary(CTML::Node("summary").AppendChild(CTML::Node("b", desc.parameters().getNameString())));
    CTML::Node mod_details("details");
    const auto desc_type = desc.type();
    if (desc_type == ParametersDescription::Type::ParametersVector) {
    } else if (desc_type == ParametersDescription::Type::Value) {
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
        } else if (subdesc_type == ParametersDescription::Type::ParametersVector) {
          item.AppendText(" vector of parameters");
          if (!subdesc.description().empty())
            item.AppendText(" defining a ").AppendChild(CTML::Node("i", subdesc.description()));
          item.AppendChild(moduleDescription(subdesc));
          const auto& vparams = desc.parameters().get<std::vector<ParametersList> >(key);
          if (!vparams.empty()) {
            CTML::Node itparams("ol");
            for (const auto& it : vparams)
              itparams.AppendChild(CTML::Node("li").AppendChild(moduleDescription(ParametersDescription(it))));
            item.AppendChild(CTML::Node("details")
                                 .AppendChild(CTML::Node("summary").AppendChild(CTML::Node("b", "Default content")))
                                 .AppendChild(CTML::Node("p").AppendChild(itparams)));
          }
        } else
          item.AppendChild(CTML::Node("i", " " + subdesc.description())).AppendChild(moduleDescription(subdesc));
        items.AppendChild(item);
      }
      if (!items.GetChildren().empty())
        mod_details.AppendChild(CTML::Node("p").AppendText("List of parameters:").AppendChild(items));
      out.AppendChild(mod_details);
    } catch (const Exception& exc) {
      exc.dump();
    }
    return out;
  }

  const std::string page_title_;
  const bool bare_, show_git_;
  CTML::Document doc_;
  CTML::Node container_;
};
REGISTER_DOCUMENTATION_GENERATOR("ctml", CTMLDocumentationGenerator);
