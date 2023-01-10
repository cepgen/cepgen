/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2023  Laurent Forthomme
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

#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Parameters.h"

namespace cepgen {
  /**
     * \brief Handler for the generic text file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
  class TextVariablesHandler : public EventExporter {
  public:
    explicit TextVariablesHandler(const ParametersList&);
    ~TextVariablesHandler();

    static ParametersDescription description();

    void initialise() override;
    void setCrossSection(double, double) override {}
    void operator<<(const Event&) override;

  private:
    std::ofstream file_;
    //--- variables definition
    const std::vector<std::string> variables_;
    const bool save_banner_, save_variables_;
    const std::string separator_;

    const utils::EventBrowser browser_;

    std::ostringstream oss_vars_;
  };

  TextVariablesHandler::TextVariablesHandler(const ParametersList& params)
      : EventExporter(params),
        file_(steer<std::string>("filename")),
        variables_(steer<std::vector<std::string> >("variables")),
        save_banner_(steer<bool>("saveBanner")),
        save_variables_(steer<bool>("saveVariables")),
        separator_(steer<std::string>("separator")) {
    //--- extract list of variables to store in output file
    oss_vars_.clear();
    std::string sep;
    for (const auto& var : variables_)
      oss_vars_ << sep << var, sep = separator_;
  }

  TextVariablesHandler::~TextVariablesHandler() {
    //--- finalisation of the output file
    file_.close();
  }

  void TextVariablesHandler::initialise() {
    if (save_banner_)
      file_ << banner("#") << "\n";
    if (save_variables_)
      file_ << "# " << oss_vars_.str() << "\n";
  }

  void TextVariablesHandler::operator<<(const Event& ev) {
    //--- write down the variables list in the file
    if (!variables_.empty()) {
      std::string sep;
      for (const auto& var : variables_)
        file_ << sep << browser_.get(ev, var), sep = separator_;
      file_ << "\n";
    }
  }

  ParametersDescription TextVariablesHandler::description() {
    auto desc = EventExporter::description();
    desc.setDescription("Text dump of variables");
    desc.add<std::string>("filename", "output.txt").setDescription("Output filename for variables dump");
    desc.add<std::vector<std::string> >("variables", {}).setDescription("List of variables to dump");
    desc.add<bool>("saveBanner", true).setDescription("Also save the boilerplate in output files?");
    desc.add<bool>("saveVariables", true).setDescription("Save the variable(s) into an output file?");
    desc.add<std::string>("separator", "\t").setDescription("Base separator in output file");
    return desc;
  }
}  // namespace cepgen

REGISTER_EXPORTER("vars", TextVariablesHandler)
