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
#include "CepGenHepMC2/HepMC2EventInterface.h"

using namespace cepgen;
using namespace HepMC;
using namespace std::string_literals;

/// Handler for the HepMC file output
/// \tparam T HepMC writer handler (format-dependent)
/// \author Laurent Forthomme <laurent.forthomme@cern.ch>
/// \date Sep 2016
template <typename T>
class HepMC2Handler : public EventExporter {
public:
  explicit HepMC2Handler(const ParametersList& params)
      : EventExporter(params), output_(new T(steer<std::string>("filename").c_str())), xs_(new GenCrossSection) {
    CG_INFO("HepMC") << "Interfacing module initialised "
                     << "for HepMC version " << HEPMC_VERSION << ".";
  }

  static ParametersDescription description() {
    auto desc = EventExporter::description();
    desc.setDescription("HepMC2 ASCII file output module");
    desc.add("filename", "output.hepmc"s).setDescription("Output filename");
    return desc;
  }

  /// Writer operator
  bool operator<<(const Event& cg_evt) override {
    CepGenEvent event(cg_evt);
    event.set_cross_section(*xs_);
    event.set_event_number(event_num_++);
    output_->write_event(&event);
    CG_DEBUG("HepMC2Handler").log([&event](auto& log) {
      log << "\n";
      event.print(log.stream());
    });
    return true;
  }
  void setCrossSection(const Value& cross_section) override {
    xs_->set_cross_section(cross_section, cross_section.uncertainty());
  }

private:
  void initialise() override {}

  const std::unique_ptr<T> output_;            ///< writer object
  const std::shared_ptr<GenCrossSection> xs_;  ///< generator cross-section and error
};

//----------------------------------------------------------------------
// Defining the various templated plugins made available by this
// specific version of HepMC (v2 and below)
//----------------------------------------------------------------------
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/IO_GenEvent.h"
typedef HepMC2Handler<IO_GenEvent> HepMC2GenEventHandler;
typedef HepMC2Handler<IO_AsciiParticles> HepMC2AsciiHandler;
REGISTER_EXPORTER("hepmc2", HepMC2GenEventHandler);
REGISTER_EXPORTER("hepmc2_ascii", HepMC2AsciiHandler);
