/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2025  Laurent Forthomme
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

#include <sstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Utils/Caller.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Value.h"
#include "CepGenPythia8/CepGenEvent.h"

using namespace std::string_literals;

namespace cepgen::pythia8 {
  /// Pythia8 handler for the LHE file output
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Sep 2016
  class LHEFEventExporter final : public EventExporter {
  public:
    explicit LHEFEventExporter(const ParametersList& params)
        : EventExporter(params),
          pythia_(new Pythia8::Pythia),
          cepgen_event_(new Pythia8::CepGenEvent),
          compress_event_(steer<bool>("compress")),
          filename_(steer<std::string>("filename")) {
      if (utils::fileExtension(filename_) == ".gz") {
#ifdef GZIP_BIN
        utils::replaceAll(filename_, ".gz", "");
#else
        CG_WARNING("pythia8:LHEFHandler")
            << "gzip compression requested, but the executable was not linked at Pythia8 wrapper compile time.";
#endif
        gzip_ = true;
      }
      if (auto file_tmp = std::ofstream(filename_); !file_tmp.is_open())
        throw CG_FATAL("pythia8:LHEFHandler") << "Failed to open output filename '" << filename_ << "' for writing.";
      cepgen_event_->openLHEF(filename_);
    }
    ~LHEFEventExporter() override {
      if (cepgen_event_)
        cepgen_event_->closeLHEF(false);  // we do not want to rewrite the init block
      if (gzip_)
#ifdef GZIP_BIN
        utils::Caller::call({GZIP_BIN, "-f", filename_});
#endif
    }

    static ParametersDescription description() {
      auto desc = EventExporter::description();
      desc.setDescription("Pythia 8-based LHEF output module");
      desc.add("compress", true);
      desc.add("filename", "output.lhe"s).setDescription("Output filename");
      return desc;
    }

    void initialise() override {
      std::ostringstream oss_init;
      oss_init << "<!--\n" << banner() << "\n-->";
      oss_init << std::endl;  // LHEF is usually not as beautifully parsed as a standard XML...
                              // we're physicists, what do you expect?
      cepgen_event_->addComments(oss_init.str());
      cepgen_event_->initialise(runParameters());
#if PYTHIA_VERSION_INTEGER < 8300
      pythia_->setLHAupPtr(cepgen_event_.get());
#else
      pythia_->setLHAupPtr(cepgen_event_);
#endif
      pythia_->settings.flag("ProcessLevel:all", false);  // we do not want Pythia to interfere...
      pythia_->settings.flag("PartonLevel:all", false);   // we do not want Pythia to interfere...
      pythia_->settings.flag("HadronLevel:all", false);   // we do not want Pythia to interfere...
      pythia_->settings.mode("Beams:frameType", 5);       // LHEF event readout
      pythia_->settings.mode("Next:numberCount", 0);      // remove some of the Pythia output
      pythia_->init();
      cepgen_event_->initLHEF();
    }

    bool operator<<(const Event& ev) override {
      cepgen_event_->feedEvent(compress_event_ ? ev : ev.compress(),
                               Pythia8::CepGenEvent::Type::centralAndFullBeamRemnants);
      pythia_->next();
      cepgen_event_->eventLHEF();
      return true;
    }
    void setCrossSection(const Value& cross_section) override {
      cepgen_event_->setCrossSection(0, cross_section, cross_section.uncertainty());
    }

  private:
    const std::unique_ptr<Pythia8::Pythia> pythia_;
    const std::shared_ptr<Pythia8::CepGenEvent> cepgen_event_;
    const bool compress_event_;
    std::string filename_;
    bool gzip_{false};
  };
}  // namespace cepgen::pythia8
using Pythia8LHEFEventExporter = cepgen::pythia8::LHEFEventExporter;
REGISTER_EXPORTER("lhef", Pythia8LHEFEventExporter);
