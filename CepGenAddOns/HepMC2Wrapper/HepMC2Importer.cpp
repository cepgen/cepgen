/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGenAddOns/HepMC2Wrapper/HepMC2EventInterface.h"

using namespace HepMC;

namespace cepgen {
  /// Handler for the HepMC file output
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Dec 2022
  class HepMC2Importer final : public EventImporter {
  public:
    /// Class constructor
    explicit HepMC2Importer(const ParametersList& params) : EventImporter(params) {
      CG_INFO("HepMC") << "Interfacing module initialised "
                       << "for HepMC version " << HEPMC_VERSION << ".";
    }

    static ParametersDescription description() {
      auto desc = EventImporter::description();
      desc.setDescription("HepMC2 ASCII file importer module");
      desc.add<std::string>("filename", "input.hepmc").setDescription("Input filename");
      return desc;
    }

    void convert(const void* in, Event& evt) const override {
      auto* evt_in = static_cast<const HepMC::GenEvent*>(in);
      evt_in->print();
      evt = Event(*static_cast<const HepMC::CepGenEvent*>(in));
    }

  private:
  };
}  // namespace cepgen
REGISTER_EVENT_IMPORTER("hepmc2", HepMC2Importer)
