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

#include <TFile.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGen/Utils/Value.h"
#include "CepGenAddOns/ROOTWrapper/ROOTTreeInfo.h"

namespace cepgen {
  /// ROOT handler for an event tree import
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Feb 2024
  class ROOTTreeImporter : public EventImporter {
  public:
    /// Class constructor
    explicit ROOTTreeImporter(const ParametersList& params)
        : EventImporter(params), file_(TFile::Open(steer<std::string>("filename").data())) {
      if (!file_)
        throw CG_FATAL("ROOTTreeImporter")
            << "Failed to load the ROOT file '" << steer<std::string>("filename") << "'.";
      run_tree_.attach(file_.get());
      evt_tree_.attach(file_.get());
    }

    static ParametersDescription description() {
      auto desc = EventImporter::description();
      desc.setDescription("ROOT TTree importer module");
      desc.add<std::string>("filename", "output.root").setDescription("Input filename");
      return desc;
    }

    bool operator>>(Event& evt) const override { return evt_tree_.next(evt); }

  private:
    void initialise() override { setCrossSection(Value{run_tree_.xsect, run_tree_.errxsect}); }

    const std::unique_ptr<TFile> file_;
    ROOT::CepGenRun run_tree_;
    mutable ROOT::CepGenEvent evt_tree_;
  };
}  // namespace cepgen
REGISTER_EVENT_IMPORTER("root_tree", ROOTTreeImporter);
