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

#include <HepMC3/GenCrossSection.h>
#include <HepMC3/Version.h>

#include <memory>

#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/Value.h"
#include "CepGenHepMC3/CepGenEvent.h"

using namespace HepMC3;
using namespace std::string_literals;

namespace cepgen::hepmc3 {
  /// Handler for the HepMC3 file output
  /// \tparam T HepMC writer handler (format-dependent)
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Sep 2016
  template <typename T>
  class EventExporter final : public cepgen::EventExporter {
  public:
    explicit EventExporter(const ParametersList& params)
        : cepgen::EventExporter(params),
          output_(new T(steer<std::string>("filename").c_str())),
          xs_(new GenCrossSection),
          run_info_(new GenRunInfo) {
      output_->set_run_info(run_info_);
      run_info_->set_weight_names({"Default"});
      CG_INFO("HepMC") << "Interfacing module initialised "
                       << "for HepMC version " << HEPMC3_VERSION << ".";
    }
    ~EventExporter() override { output_->close(); }

    static ParametersDescription description() {
      auto desc = cepgen::EventExporter::description();
      desc.setDescription("HepMC3 ASCII file output module");
      desc.add("filename", "output.hepmc"s).setDescription("Output filename");
      return desc;
    }

    bool operator<<(const Event& cg_event) override {
      CepGenEvent event(cg_event);
      event.set_cross_section(xs_);
      event.set_run_info(run_info_);
      event.set_event_number(event_num_++);
      output_->write_event(event);
      return !output_->failed();
    }
    void setCrossSection(const Value& cross_section) override {
      xs_->set_cross_section(cross_section, cross_section.uncertainty());
    }

  private:
    void initialise() override {}

    const std::unique_ptr<T> output_;             ///< writer object
    const std::shared_ptr<GenCrossSection> xs_;   ///< generator cross-section and error
    const std::shared_ptr<GenRunInfo> run_info_;  ///< auxiliary information on run
  };
}  // namespace cepgen::hepmc3
using HepMC3Handler = cepgen::hepmc3::EventExporter;
//----------------------------------------------------------------------
// Defining the various templated plugins made available by this
// specific version of HepMC
//----------------------------------------------------------------------
#include <HepMC3/WriterAscii.h>
#include <HepMC3/WriterHEPEVT.h>
using HepMC3AsciiHandler = HepMC3Handler<WriterAscii>;
using HepMC3HEPEVTHandler = HepMC3Handler<WriterHEPEVT>;
REGISTER_EXPORTER("hepmc", HepMC3AsciiHandler);
REGISTER_EXPORTER("hepevt", HepMC3HEPEVTHandler);

#if HEPMC3_VERSION_CODE >= 3001000
#include <HepMC3/WriterAsciiHepMC2.h>
using HepMC3HepMC2Handler = HepMC3Handler<WriterAsciiHepMC2>;
REGISTER_EXPORTER("hepmc3_hepmc2", HepMC3HepMC2Handler);
#if HEPMC3_VERSION_CODE >= 3002005 && HEPMC3_USE_COMPRESSION
#include <HepMC3/WriterGZ.h>
using HepMC3AsciiZHandler = HepMC3Handler<WriterGZ<WriterAscii> >;
using HepMC3HEPEVTZHandler = HepMC3Handler<WriterGZ<WriterHEPEVT> >;
using HepMC3AsciiLZMAHandler = HepMC3Handler<WriterGZ<WriterAscii, Compression::lzma> >;
using HepMC3HEPEVTLZMAHandler = HepMC3Handler<WriterGZ<WriterHEPEVT, Compression::lzma> >;
using HepMC3AsciiBZ2Handler = HepMC3Handler<WriterGZ<WriterAscii, Compression::bz2> >;
using HepMC3HEPEVTBZ2Handler = HepMC3Handler<WriterGZ<WriterHEPEVT, Compression::bz2> >;
REGISTER_EXPORTER("hepmc_z", HepMC3AsciiZHandler);
REGISTER_EXPORTER("hepevt_z", HepMC3HEPEVTZHandler);
REGISTER_EXPORTER("hepmc_lzma", HepMC3AsciiLZMAHandler);
REGISTER_EXPORTER("hepevt_lzma", HepMC3HEPEVTLZMAHandler);
REGISTER_EXPORTER("hepmc_bz2", HepMC3AsciiBZ2Handler);
REGISTER_EXPORTER("hepevt_bz2", HepMC3HEPEVTBZ2Handler);
#endif
#endif

#ifdef HEPMC3_ROOTIO
#include <HepMC3/WriterRoot.h>
#include <HepMC3/WriterRootTree.h>
using HepMC3RootHandler = HepMC3Handler<WriterRoot>;
using HepMC3RootTreeHandler = HepMC3Handler<WriterRootTree>;
REGISTER_EXPORTER("hepmc_root", HepMC3RootHandler);
REGISTER_EXPORTER("hepmc_root_tree", HepMC3RootTreeHandler);
#endif

#ifdef HEPMC3_EXTRA_PLUGINS
#include <ConvertExample/include/WriterDOT.h>
#include <ConvertExample/include/WriterRootTreeOPAL.h>
using HepMC3DOTHandler = HepMC3Handler<WriterDOT>;
using HepMC3RootTreeOPALHandler = HepMC3Handler<WriterRootTreeOPAL>;
REGISTER_EXPORTER("hepmc_dot", HepMC3DOTHandler);
REGISTER_EXPORTER("hepmc_root_tree_opal", HepMC3RootTreeOPALHandler);
#endif
