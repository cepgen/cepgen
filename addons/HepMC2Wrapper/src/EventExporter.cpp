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

#include <HepMC/Version.h>

#include <memory>

#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/Value.h"
#include "CepGenHepMC2/CepGenEvent.h"

using namespace HepMC;
using namespace std::string_literals;

namespace cepgen::hepmc2 {
  /// Handler for the HepMC file output
  /// \tparam T HepMC writer handler (format-dependent)
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Sep 2016
  template <typename T>
  class EventExporter final : public cepgen::EventExporter {
  public:
    explicit EventExporter(const ParametersList& params)
        : cepgen::EventExporter(params),
          output_(new T(steer<std::string>("filename").c_str())),
          cross_section_(new GenCrossSection) {
      CG_INFO("HepMC") << "Interfacing module initialised "
                       << "for HepMC version " << HEPMC_VERSION << ".";
    }

    static ParametersDescription description() {
      auto desc = cepgen::EventExporter::description();
      desc.setDescription("HepMC2 ASCII file output module");
      desc.add("filename", "output.hepmc"s).setDescription("Output filename");
      return desc;
    }

    /// Writer operator
    bool operator<<(const Event& cepgen_event) override {
      CepGenEvent event(cepgen_event);
      event.set_cross_section(*cross_section_);
      event.set_event_number(event_num_++);
      output_->write_event(&event);
      CG_DEBUG("HepMC2Handler").log([&event](auto& log) {
        log << "\n";
        event.print(log.stream());
      });
      return true;
    }
    void setCrossSection(const Value& cross_section) override {
      cross_section_->set_cross_section(cross_section, cross_section.uncertainty());
    }

  private:
    void initialise() override {}

    const std::unique_ptr<T> output_;                       ///< writer object
    const std::shared_ptr<GenCrossSection> cross_section_;  ///< generator cross-section and error
  };
}  // namespace cepgen::hepmc2
using cepgen::hepmc2::EventExporter;
//----------------------------------------------------------------------
// Defining the various templated plugins made available by this
// specific version of HepMC (v2 and below)
//----------------------------------------------------------------------
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/IO_GenEvent.h"
using HepMC2GenEventHandler = EventExporter<IO_GenEvent>;
using HepMC2AsciiHandler = EventExporter<IO_AsciiParticles>;
REGISTER_EXPORTER("hepmc2", HepMC2GenEventHandler);
REGISTER_EXPORTER("hepmc2_ascii", HepMC2AsciiHandler);
