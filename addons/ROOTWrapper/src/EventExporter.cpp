/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2014-2025  Laurent Forthomme
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
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Value.h"
#include "CepGen/Version.h"
#include "CepGenRoot/ROOTTreeInfo.h"

using namespace std::string_literals;

namespace cepgen::root {
  /// Handler for the storage of events in a ROOT format
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date 27 Jan 2014
  class EventExporter final : public cepgen::EventExporter {
  public:
    explicit EventExporter(const ParametersList& params)
        : cepgen::EventExporter(params),
          filename_(steer<std::string>("filename")),
          compress_(steer<bool>("compress")),
          file_(TFile::Open(filename_.data(), "recreate")) {
      if (!file_->IsOpen())
        throw CG_FATAL("root:EventExporter") << "Failed to create the output file!";
    }
    ~EventExporter() override {
      run_tree_.fill();
      file_->Write();
    }

    static ParametersDescription description() {
      auto desc = cepgen::EventExporter::description();
      desc.setDescription("ROOT TTree storage module");
      desc.add("filename", "output.root"s).setDescription("Output filename");
      desc.add("compress", false).setDescription("Compress the event content? (merge down two-parton system)");
      desc.add("autoFilename", false).setDescription("automatically generate the output filename");
      return desc;
    }

    bool operator<<(const Event& event) override {
      event_tree_.fill(event, compress_);
      run_tree_.num_events += 1;
      return true;
    }
    void setCrossSection(const Value& cross_section) override {
      run_tree_.xsect = cross_section;
      run_tree_.errxsect = cross_section.uncertainty();
    }

  private:
    void initialise() override {
      if (steer<bool>("autoFilename")) {
        auto filename = generateFilename();
        CG_INFO("root:EventExporter") << "Output ROOT filename automatically set to '" << filename << "'.";
        if (file_.reset(TFile::Open(filename.data(), "recreate")); !file_->IsOpen())
          throw CG_FATAL("root:EventExporter") << "Failed to create the output file!";
      }
      run_tree_.create();
      event_tree_.create();
      run_tree_.litigious_events = 0;
      if (runParameters().hasProcess()) {
        run_tree_.sqrt_s = runParameters().kinematics().incomingBeams().sqrtS();
        run_tree_.process_name = runParameters().processName();
        run_tree_.process_parameters = runParameters().process().parameters().serialise();
      }
    }
    std::string generateFilename() const {
      std::string event_modifiers, process_mode;
      for (const auto& mod : runParameters().eventModifiersSequence())
        event_modifiers += (event_modifiers.empty() ? "" : "-") + mod->name();
      const auto symmetrise = runParameters().process().parameters().get<bool>("symmetrise");
      const auto sf_info =
          utils::sanitise(runParameters().process().kinematics().incomingBeams().structureFunctions().serialise());
      switch (runParameters().process().kinematics().incomingBeams().mode()) {
        case mode::Kinematics::ElasticElastic:
          process_mode = "el";
          break;
        case mode::Kinematics::InelasticElastic:
          process_mode = symmetrise ? "sd" : "sdie_" + sf_info;
          break;
        case mode::Kinematics::ElasticInelastic:
          process_mode = symmetrise ? "sd" : "sdei_" + sf_info;
          break;
        case mode::Kinematics::InelasticInelastic:
          process_mode = "dd_" + sf_info;
          break;
        case mode::Kinematics::invalid:
          break;
      }
      return utils::format("cepgen%s_%s_%s_%gTeV%s.root",
                           utils::sanitise(version::tag).data(),
                           runParameters().processName().data(),
                           process_mode.data(),
                           runParameters().kinematics().incomingBeams().sqrtS() / 1000.,
                           event_modifiers.data());
    }

    const std::string filename_;
    const bool compress_;
    std::unique_ptr<TFile> file_;
    ROOT::CepGenRun run_tree_;
    ROOT::CepGenEvent event_tree_;
  };
}  // namespace cepgen::root
using ROOTEventExporter = cepgen::root::EventExporter;
REGISTER_EXPORTER("root_tree", ROOTEventExporter);
