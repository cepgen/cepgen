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

#include <Rivet/Analysis.hh>
#include <Rivet/AnalysisHandler.hh>
#include <Rivet/Tools/RivetPaths.hh>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/HepMC2Wrapper/HepMC2EventInterface.h"

namespace cepgen {
  /**
     * \brief Handler for the Rivet analysis framework
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Aug 2019
     */
  class RivetAnalysisHandler final : public EventExporter {
  public:
    explicit RivetAnalysisHandler(const ParametersList&);
    ~RivetAnalysisHandler();

    static ParametersDescription description() {
      auto desc = EventExporter::description();
      desc.setDescription("Rivet analysis handler");
      desc.add<std::string>("filename", "output.rivet.yoda");
      desc.add<std::vector<std::string> >("analyses", {});
      return desc;
    }

    void initialise() override;
    void setCrossSection(double xsec, double err_xsec) override { xsec_ = xsec, err_xsec_ = err_xsec; }
    void operator<<(const Event&) override;

  private:
    void* enginePtr() override { return (void*)rivet_.get(); }

    std::unique_ptr<Rivet::AnalysisHandler> rivet_;
    const std::string filename_;
    const std::vector<std::string> analyses_;

    double xsec_{0.}, err_xsec_{0.};
  };

  RivetAnalysisHandler::RivetAnalysisHandler(const ParametersList& params)
      : EventExporter(params),
        rivet_(new Rivet::AnalysisHandler("CepGen")),
        filename_(steer<std::string>("filename")),
        analyses_(steer<std::vector<std::string> >("analyses")) {
    if (analyses_.empty())
      throw CG_FATAL("RivetAnalysisHandler") << "At least one analysis is required!";
    for (const auto& path : params.get<std::vector<std::string> >("paths"))
      Rivet::addAnalysisLibPath(path);
    rivet_->addAnalyses(analyses_);
    if (analyses_.size() != rivet_->analysesMap().size())
      throw CG_FATAL("RivetAnalysisHandler") << "Rivet failed to find all analyses requested!\n\t"
                                             << "You may used `rivet --list-analyses` to dump a full list.";
  }

  RivetAnalysisHandler::~RivetAnalysisHandler() {
    rivet_->setCrossSection({xsec_, err_xsec_});
    rivet_->finalize();
    rivet_->writeData(filename_);
  }

  void RivetAnalysisHandler::initialise() {
    if (!runParameters().hasProcess())
      throw CG_FATAL("RivetAnalysisHandler") << "No process defined!";
    if (!runParameters().process().hasEvent())
      throw CG_FATAL("RivetAnalysisHandler")
          << "Process \"" << runParameters().processName() << "\" has no event content!";
    rivet_->init(HepMC::CepGenEvent(runParameters().process().event()));
  }

  void RivetAnalysisHandler::operator<<(const Event& ev) {
    HepMC::CepGenEvent hepmc_evt(ev);
    try {
      rivet_->analyze(hepmc_evt);
    } catch (const YODA::Exception& err) {
      CG_WARNING("RivetAnalysisHandler") << "Rivet/YODA encountered the following exception:\n\t" << err.what();
    }
  }
}  // namespace cepgen
REGISTER_EXPORTER("rivet", RivetAnalysisHandler);
