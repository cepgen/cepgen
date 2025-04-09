/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGenHepMC3/CepGenEvent.h"

using namespace std::string_literals;

namespace cepgen::hepmc3 {
  /// Handler for the HepMC file output
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Feb 2024
  template <typename T>
  class EventImporter final : public cepgen::EventImporter {
  public:
    /// Class constructor
    explicit EventImporter(const ParametersList& params)
        : cepgen::EventImporter(params), reader_(new T(steer<std::string>("filename"))) {
      if (!reader_)
        throw CG_FATAL("hepmc3:EventImporter") << "Failed to initialise HepMC reader.";
      CG_INFO("hepmc3:EventImporter") << "Interfacing module initialised "
                                      << "for HepMC version " << HEPMC3_VERSION << " and HepMC ASCII file '"
                                      << steer<std::string>("filename") << "'.";
    }

    bool operator>>(Event& evt) override {
      HepMC3::GenEvent event;
      if (!reader_->read_event(event))
        return false;
      if (!cross_section_retrieved_) {
        if (const auto xsec = event.cross_section(); xsec)
          setCrossSection(Value{xsec->xsec(), xsec->xsec_err()});
        cross_section_retrieved_ = true;
      }
      CG_DEBUG("hepmc:EventImporter").log([&event](auto& log) { HepMC3::Print::content(log.stream(), event); });
      evt = Event(static_cast<const HepMC3::CepGenEvent&>(event));
      return true;
    }

    static ParametersDescription description() {
      auto desc = cepgen::EventImporter::description();
      desc.setDescription("HepMC3 ASCII file importer module");
      desc.add("filename", "input.hepmc"s).setDescription("Input filename");
      return desc;
    }

  private:
    void initialise() override {}
    const std::unique_ptr<T> reader_;
    bool cross_section_retrieved_{false};
  };
}  // namespace cepgen::hepmc3
template <typename T>
using HepMC3EventImporter = cepgen::hepmc3::EventImporter<T>;
using HepMC3ImporterASCII = HepMC3EventImporter<HepMC3::ReaderAscii>;
using HepMC3ImporterHEPEVT = HepMC3EventImporter<HepMC3::ReaderHEPEVT>;
REGISTER_EVENT_IMPORTER("hepmc", HepMC3ImporterASCII);
//REGISTER_EVENT_IMPORTER("hepevt", HepMC3ImporterHEPEVT); // HEPEVT input is still very shaky, disabling it by default
#if HEPMC3_VERSION_CODE >= 3001000
#include <HepMC3/ReaderAsciiHepMC2.h>
using HepMC3ImporterHepMC2 = HepMC3EventImporter<HepMC3::ReaderAsciiHepMC2>;
REGISTER_EVENT_IMPORTER("hepmc3_hepmc2", HepMC3ImporterHepMC2);
/*#if HEPMC3_VERSION_CODE >= 3002005 && HEPMC3_USE_COMPRESSION
#include <HepMC3/ReaderGZ.h>
using HepMC3ImporterAsciiZ = HepMC3EventImporter<HepMC3::ReaderGZ<HepMC3::ReaderAscii> >;
using HepMC3ImporterHEPEVTZ = HepMC3EventImporter<HepMC3::ReaderGZ<HepMC3::ReaderHEPEVT> >;
REGISTER_EVENT_IMPORTER("hepmc_z", HepMC3ImporterAsciiZ);
REGISTER_EVENT_IMPORTER("hepevt_z", HepMC3ImporterHEPEVTZ);
#endif*/
#endif

/*#ifdef HEPMC3_ROOTIO
#include <HepMC3/ReaderRoot.h>
#include <HepMC3/ReaderRootTree.h>
using HepMC3ImporterRoot = HepMC3EventImporter<HepMC3::ReaderRoot>;
using HepMC3ImporterRootTree = HepMC3EventImporter<HepMC3::ReaderRootTree>;
REGISTER_EVENT_IMPORTER("hepmc_root", HepMC3ImporterRoot);
REGISTER_EVENT_IMPORTER("hepmc_root_tree", HepMC3ImporterRootTree);
#endif*/

#ifdef HEPMC3_EXTRA_PLUGINS
#include <ConvertExample/include/ReaderDOT.h>
#include <ConvertExample/include/ReaderRootTreeOPAL.h>
using HepMC3ImporterDOT = HepMC3EventImporter<HepMC3::ReaderDOT>;
using HepMC3ImporterRootTreeOPAL = HepMC3EventImporter<HepMC3::ReaderRootTreeOPAL>;
REGISTER_EVENT_IMPORTER("hepmc_dot", HepMC3ImporterDOT);
REGISTER_EVENT_IMPORTER("hepmc_root_tree_opal", HepMC3ImporterRootTreeOPAL);
#endif
