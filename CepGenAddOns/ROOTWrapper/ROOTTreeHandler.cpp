/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#include <sstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Parameters.h"
#include "CepGenAddOns/ROOTWrapper/ROOTTreeInfo.h"

namespace cepgen {
  namespace io {
    /**
     * Handler for the storage of events in a ROOT format
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date 27 Jan 2014
     */
    class ROOTTreeHandler : public ExportModule {
    public:
      /// Class constructor
      explicit ROOTTreeHandler(const ParametersList&);
      ~ROOTTreeHandler();

      static ParametersDescription description();

      void initialise(const Parameters&) override;
      /// Writer operator
      void operator<<(const Event&) override;
      void setCrossSection(double, double) override;

    private:
      TFile file_;
      const bool compress_;
      ROOT::CepGenRun run_tree_;
      ROOT::CepGenEvent evt_tree_;
    };

    ROOTTreeHandler::ROOTTreeHandler(const ParametersList& params)
        : ExportModule(params),
          file_(steer<std::string>("filename").c_str(), "recreate"),
          compress_(steer<bool>("compress")) {
      if (!file_.IsOpen())
        throw CG_FATAL("ROOTTreeHandler") << "Failed to create the output file!";
      run_tree_.create();
      evt_tree_.create();
    }

    ROOTTreeHandler::~ROOTTreeHandler() {
      run_tree_.fill();
      file_.Write();
    }

    void ROOTTreeHandler::initialise(const Parameters& params) {
      run_tree_.litigious_events = 0;
      run_tree_.sqrt_s = params.kinematics().incomingBeams().sqrtS();
      run_tree_.process_name = params.processName();
    }

    void ROOTTreeHandler::operator<<(const Event& ev) {
      evt_tree_.fill(ev, compress_);
      run_tree_.num_events += 1;
    }

    void ROOTTreeHandler::setCrossSection(double cross_section, double cross_section_err) {
      run_tree_.xsect = cross_section;
      run_tree_.errxsect = cross_section_err;
    }

    ParametersDescription ROOTTreeHandler::description() {
      auto desc = ExportModule::description();
      desc.setDescription("ROOT TTree storage module");
      desc.add<std::string>("filename", "output.root").setDescription("Output filename");
      desc.add<bool>("compress", false).setDescription("Compress the event content? (merge down two-parton system)");
      return desc;
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("root_tree", ROOTTreeHandler)
