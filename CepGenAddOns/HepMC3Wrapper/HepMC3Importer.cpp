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

#include <HepMC3/GenEvent.h>
#include <HepMC3/Print.h>
#include <HepMC3/ReaderAscii.h>
#include <HepMC3/ReaderHEPEVT.h>
#include <HepMC3/Version.h>

#include <memory>

#include "CepGen/Core/Exception.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGenAddOns/HepMC3Wrapper/HepMC3EventInterface.h"

namespace cepgen {
  /// Handler for the HepMC file output
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Feb 2024
  template <typename T>
  class HepMC3Importer final : public EventImporter {
  public:
    /// Class constructor
    explicit HepMC3Importer(const ParametersList& params)
        : EventImporter(params), reader_(new T(steer<std::string>("filename"))) {
      if (!reader_)
        throw CG_FATAL("HepMC3Importer") << "Failed to initialise HepMC reader.";
      CG_INFO("HepMC3Importer") << "Interfacing module initialised "
                                << "for HepMC version " << HEPMC3_VERSION << " and HepMC ASCII file '"
                                << steer<std::string>("filename") << "'.";
    }

    bool operator>>(Event& evt) const override {
      HepMC3::GenEvent event;
      if (!reader_->read_event(event))
        return false;
      CG_DEBUG("HepMC3Importer").log([&event](auto& log) { HepMC3::Print::content(log.stream(), event); });
      evt = Event(static_cast<const HepMC3::CepGenEvent&>(event));
      return true;
    }

    static ParametersDescription description() {
      auto desc = EventImporter::description();
      desc.setDescription("HepMC3 ASCII file importer module");
      desc.add<std::string>("filename", "input.hepmc").setDescription("Input filename");
      return desc;
    }

  private:
    void initialise() override {}
    const std::unique_ptr<T> reader_;
  };
}  // namespace cepgen
typedef cepgen::HepMC3Importer<HepMC3::ReaderAscii> HepMC3ImporterASCII;
typedef cepgen::HepMC3Importer<HepMC3::ReaderHEPEVT> HepMC3ImporterHEPEVT;
REGISTER_EVENT_IMPORTER("hepmc", HepMC3ImporterASCII);
REGISTER_EVENT_IMPORTER("hepevt", HepMC3ImporterHEPEVT);
