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

#include <TFile.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGen/Utils/Value.h"
#include "CepGenRoot/ROOTTreeInfo.h"

using namespace std::string_literals;

namespace cepgen::root {
  /// ROOT handler for an event tree import
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Feb 2024
  class EventImporter : public cepgen::EventImporter {
  public:
    explicit EventImporter(const ParametersList& params)
        : cepgen::EventImporter(params), file_(TFile::Open(steer<std::string>("filename").data())) {
      if (!file_)
        throw CG_FATAL("root::EventImporter")
            << "Failed to load the ROOT file '" << steer<std::string>("filename") << "'.";
      run_tree_.attach(file_.get());
      event_tree_.attach(file_.get());
    }

    static ParametersDescription description() {
      auto desc = cepgen::EventImporter::description();
      desc.setDescription("ROOT TTree importer module");
      desc.add("filename", "output.root"s).setDescription("Input filename");
      return desc;
    }

    bool operator>>(Event& event) override { return event_tree_.next(event); }

  private:
    void initialise() override { setCrossSection(Value{run_tree_.xsect, run_tree_.errxsect}); }

    const std::unique_ptr<TFile> file_;
    ROOT::CepGenRun run_tree_;
    ROOT::CepGenEvent event_tree_;
  };
}  // namespace cepgen::root
using ROOTEventImporter = cepgen::root::EventImporter;
REGISTER_EVENT_IMPORTER("root_tree", ROOTEventImporter);
