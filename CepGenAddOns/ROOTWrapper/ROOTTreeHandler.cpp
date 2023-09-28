/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Value.h"
#include "CepGen/Version.h"
#include "CepGenAddOns/ROOTWrapper/ROOTTreeInfo.h"

namespace cepgen {
  /**
     * Handler for the storage of events in a ROOT format
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date 27 Jan 2014
     */
  class ROOTTreeHandler : public EventExporter {
  public:
    /// Class constructor
    explicit ROOTTreeHandler(const ParametersList&);
    ~ROOTTreeHandler();

    static ParametersDescription description();

    void initialise() override;
    /// Writer operator
    void operator<<(const Event&) override;
    void setCrossSection(const Value&) override;

  private:
    std::string generateFilename() const;

    const std::string filename_;
    const bool compress_;
    const bool auto_filename_;
    std::unique_ptr<TFile> file_;
    ROOT::CepGenRun run_tree_;
    ROOT::CepGenEvent evt_tree_;
  };

  ROOTTreeHandler::ROOTTreeHandler(const ParametersList& params)
      : EventExporter(params),
        filename_(steer<std::string>("filename")),
        compress_(steer<bool>("compress")),
        auto_filename_(steer<bool>("autoFilename")) {}

  ROOTTreeHandler::~ROOTTreeHandler() {
    run_tree_.fill();
    file_->Write();
  }

  void ROOTTreeHandler::initialise() {
    auto filename = filename_;
    if (auto_filename_) {
      filename = generateFilename();
      CG_INFO("ROOTTreeHandler") << "Output ROOT filename automatically set to '" << filename << "'.";
    }
    file_.reset(TFile::Open(filename.data(), "recreate"));
    if (!file_->IsOpen())
      throw CG_FATAL("ROOTTreeHandler") << "Failed to create the output file!";
    run_tree_.create();
    evt_tree_.create();
    run_tree_.litigious_events = 0;
    run_tree_.sqrt_s = runParameters().kinematics().incomingBeams().sqrtS();
    run_tree_.process_name = runParameters().processName();
    run_tree_.process_parameters = runParameters().process().parameters().serialise();
  }

  void ROOTTreeHandler::operator<<(const Event& ev) {
    evt_tree_.fill(ev, compress_);
    run_tree_.num_events += 1;
  }

  void ROOTTreeHandler::setCrossSection(const Value& cross_section) {
    run_tree_.xsect = cross_section;
    run_tree_.errxsect = cross_section.uncertainty();
  }

  std::string ROOTTreeHandler::generateFilename() const {
    std::string evt_mods, proc_mode;
    for (const auto& mod : runParameters().eventModifiersSequence())
      evt_mods += (evt_mods.empty() ? "" : "-") + mod->name();
    const auto symm = runParameters().process().parameters().get<bool>("symmetrise");
    const auto sf_info =
        utils::sanitise(runParameters().process().kinematics().incomingBeams().structureFunctions().serialise());
    switch (runParameters().process().kinematics().incomingBeams().mode()) {
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
                         runParameters().processName().data(),
                         proc_mode.data(),
                         runParameters().kinematics().incomingBeams().sqrtS() / 1000.,
                         evt_mods.data());
  }

  ParametersDescription ROOTTreeHandler::description() {
    auto desc = EventExporter::description();
    desc.setDescription("ROOT TTree storage module");
    desc.add<std::string>("filename", "output.root").setDescription("Output filename");
    desc.add<bool>("compress", false).setDescription("Compress the event content? (merge down two-parton system)");
    desc.add<bool>("autoFilename", false).setDescription("automatically generate the output filename");
    return desc;
  }
}  // namespace cepgen

REGISTER_EXPORTER("root_tree", ROOTTreeHandler);
