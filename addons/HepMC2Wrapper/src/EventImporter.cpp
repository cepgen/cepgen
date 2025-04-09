/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#include <HepMC/GenEvent.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/Version.h>

#include <memory>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGenHepMC2/CepGenEvent.h"

using namespace std::string_literals;

namespace cepgen::hepmc2 {
  /// Handler for the HepMC file output
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Dec 2022
  class EventImporter final : public cepgen::EventImporter {
  public:
    /// Class constructor
    explicit EventImporter(const ParametersList& params)
        : cepgen::EventImporter(params), reader_(new HepMC::IO_GenEvent(steer<std::string>("filename"), std::ios::in)) {
      if (!reader_)
        throw CG_FATAL("HepMC2Importer") << "Failed to initialise HepMCv2 reader.";
      CG_INFO("HepMC2Importer") << "Interfacing module initialised "
                                << "for HepMC version " << HEPMC_VERSION << " and HepMC ASCII file '"
                                << steer<std::string>("filename") << "' with I/O state " << reader_->rdstate() << ".";
    }

    bool operator>>(Event& evt) override {
      HepMC::GenEvent event;
      if (!reader_->fill_next_event(&event))
        return false;
      if (!cross_section_retrieved_) {
        if (const auto xsec = event.cross_section(); xsec)
          setCrossSection(Value{xsec->cross_section(), xsec->cross_section_error()});
        cross_section_retrieved_ = true;
      }
      CG_DEBUG("HepMC2Importer").log([&event](auto& log) { event.print(log.stream()); });
      evt = Event(static_cast<const HepMC::CepGenEvent&>(event));
      return true;
    }

    static ParametersDescription description() {
      auto desc = cepgen::EventImporter::description();
      desc.setDescription("HepMC2 ASCII file importer module");
      desc.add("filename", "input.hepmc"s).setDescription("Input filename");
      return desc;
    }

  private:
    void initialise() override {}
    const std::unique_ptr<HepMC::IO_GenEvent> reader_;
    bool cross_section_retrieved_{false};
  };
}  // namespace cepgen::hepmc2
using cepgen::hepmc2::EventImporter;
REGISTER_EVENT_IMPORTER("hepmc2", EventImporter);
