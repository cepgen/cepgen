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

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"
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
      std::string generateFilename(const Parameters&) const;

      const std::string filename_;
      const bool compress_;
      const bool auto_filename_;
      std::unique_ptr<TFile> file_;
      ROOT::CepGenRun run_tree_;
      ROOT::CepGenEvent evt_tree_;
    };

    ROOTTreeHandler::ROOTTreeHandler(const ParametersList& params)
        : ExportModule(params),
          filename_(steer<std::string>("filename")),
          compress_(steer<bool>("compress")),
          auto_filename_(steer<bool>("autoFilename")) {}

    ROOTTreeHandler::~ROOTTreeHandler() {
      run_tree_.fill();
      file_->Write();
    }

    void ROOTTreeHandler::initialise(const Parameters& params) {
      auto filename = filename_;
      if (auto_filename_) {
        filename = generateFilename(params);
        CG_INFO("ROOTTreeHandler") << "Output ROOT filename automatically set to '" << filename << "'.";
      }
      file_.reset(TFile::Open(filename.data(), "recreate"));
      if (!file_->IsOpen())
        throw CG_FATAL("ROOTTreeHandler") << "Failed to create the output file!";
      run_tree_.create();
      evt_tree_.create();
      run_tree_.litigious_events = 0;
      run_tree_.sqrt_s = params.kinematics().incomingBeams().sqrtS();
      run_tree_.process_name = params.processName();
      run_tree_.process_parameters = params.process().parameters().serialise();
    }

    void ROOTTreeHandler::operator<<(const Event& ev) {
      evt_tree_.fill(ev, compress_);
      run_tree_.num_events += 1;
    }

    void ROOTTreeHandler::setCrossSection(double cross_section, double cross_section_err) {
      run_tree_.xsect = cross_section;
      run_tree_.errxsect = cross_section_err;
    }

    std::string ROOTTreeHandler::generateFilename(const Parameters& params) const {
      std::string evt_mods, proc_mode;
      for (const auto& mod : params.eventModifiersSequence())
        evt_mods += (evt_mods.empty() ? "" : "-") + mod->name();
      const auto symm = params.process().parameters().get<bool>("symmetrise");
      const auto sf_info = utils::sanitise(strfun::StructureFunctionsFactory::get().describe(
          params.process().kinematics().incomingBeams().structureFunctions()->name()));
      switch (params.process().kinematics().incomingBeams().mode()) {
        case mode::Kinematics::ElasticElastic:
          proc_mode = "el";
          break;
        case mode::Kinematics::InelasticElastic:
          proc_mode = symm ? "sd" : "sdie_" + sf_info;
          break;
        case mode::Kinematics::ElasticInelastic:
          proc_mode = symm ? "sd" : "sdei_" + sf_info;
          break;
        case mode::Kinematics::InelasticInelastic:
          proc_mode = "dd_" + sf_info;
          break;
        case mode::Kinematics::invalid:
          break;
      }
      return utils::format("cepgen%s_%s_%s_%gTeV%s.root",
                           utils::sanitise(version::tag).data(),
                           params.processName().data(),
                           proc_mode.data(),
                           params.kinematics().incomingBeams().sqrtS() / 1000.,
                           evt_mods.data());
    }

    ParametersDescription ROOTTreeHandler::description() {
      auto desc = ExportModule::description();
      desc.setDescription("ROOT TTree storage module");
      desc.add<std::string>("filename", "output.root").setDescription("Output filename");
      desc.add<bool>("compress", false).setDescription("Compress the event content? (merge down two-parton system)");
      desc.add<bool>("autoFilename", false).setDescription("automatically generate the output filename");
      return desc;
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("root_tree", ROOTTreeHandler)
