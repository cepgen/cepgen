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

#include <HepMC3/GenCrossSection.h>
#include <HepMC3/Version.h>

#include <memory>

#include "CepGen/Core/EventExporter.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Parameters.h"
#include "CepGenAddOns/HepMC3Wrapper/HepMC3EventInterface.h"

using namespace HepMC3;

namespace cepgen {
  /// Handler for the HepMC3 file output
  /// \tparam T HepMC writer handler (format-dependent)
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Sep 2016
  template <typename T>
  class HepMC3Handler : public EventExporter {
  public:
    /// Class constructor
    explicit HepMC3Handler(const ParametersList&);
    ~HepMC3Handler();

    static ParametersDescription description();

    void initialise() override {}
    /// Writer operator
    void operator<<(const Event&) override;
    void setCrossSection(double, double) override;

  private:
    /// Writer object
    std::unique_ptr<T> output_;
    /// Generator cross section and error
    std::shared_ptr<GenCrossSection> xs_;
    /// Auxiliary information on run
    std::shared_ptr<GenRunInfo> run_info_;
  };

  template <typename T>
  HepMC3Handler<T>::HepMC3Handler(const ParametersList& params)
      : EventExporter(params),
        output_(new T(steer<std::string>("filename").c_str())),
        xs_(new GenCrossSection),
        run_info_(new GenRunInfo) {
    output_->set_run_info(run_info_);
    run_info_->set_weight_names({"Default"});
    CG_INFO("HepMC") << "Interfacing module initialised "
                     << "for HepMC version " << HEPMC3_VERSION << ".";
  }

  template <typename T>
  HepMC3Handler<T>::~HepMC3Handler() {
    output_->close();
  }

  template <typename T>
  void HepMC3Handler<T>::operator<<(const Event& evt) {
    CepGenEvent event(evt);
    // general information
    event.set_cross_section(xs_);
    event.set_run_info(run_info_);
    event.set_event_number(event_num_++);
    output_->write_event(event);
  }

  template <typename T>
  void HepMC3Handler<T>::setCrossSection(double cross_section, double cross_section_err) {
    xs_->set_cross_section(cross_section, cross_section_err);
  }

  template <typename T>
  ParametersDescription HepMC3Handler<T>::description() {
    auto desc = EventExporter::description();
    desc.setDescription("HepMC3 ASCII file output module");
    desc.add<std::string>("filename", "output.hepmc").setDescription("Output filename");
    return desc;
  }
}  // namespace cepgen

//----------------------------------------------------------------------
// Defining the various templated plugins made available by this
// specific version of HepMC
//----------------------------------------------------------------------

#include <HepMC3/WriterAscii.h>
#include <HepMC3/WriterHEPEVT.h>
typedef cepgen::HepMC3Handler<WriterAscii> HepMC3AsciiHandler;
typedef cepgen::HepMC3Handler<WriterHEPEVT> HepMC3HEPEVTHandler;
REGISTER_EXPORTER("hepmc", HepMC3AsciiHandler)
REGISTER_EXPORTER("hepevt", HepMC3HEPEVTHandler)

#if HEPMC3_VERSION_CODE >= 3001000
#include <HepMC3/WriterAsciiHepMC2.h>
typedef cepgen::HepMC3Handler<WriterAsciiHepMC2> HepMC3HepMC2Handler;
REGISTER_EXPORTER("hepmc3_hepmc2", HepMC3HepMC2Handler)
#if HEPMC3_VERSION_CODE >= 3002005 && HEPMC3_USE_COMPRESSION
#include <HepMC3/WriterGZ.h>
typedef cepgen::HepMC3Handler<WriterGZ<WriterAscii, Compression::z> > HepMC3AsciiZHandler;
typedef cepgen::HepMC3Handler<WriterGZ<WriterHEPEVT, Compression::z> > HepMC3HEPEVTZHandler;
typedef cepgen::HepMC3Handler<WriterGZ<WriterAscii, Compression::lzma> > HepMC3AsciiLZMAHandler;
typedef cepgen::HepMC3Handler<WriterGZ<WriterHEPEVT, Compression::lzma> > HepMC3HEPEVTLZMAHandler;
typedef cepgen::HepMC3Handler<WriterGZ<WriterAscii, Compression::bz2> > HepMC3AsciiBZ2Handler;
typedef cepgen::HepMC3Handler<WriterGZ<WriterHEPEVT, Compression::bz2> > HepMC3HEPEVTBZ2Handler;
REGISTER_EXPORTER("hepmc_z", HepMC3AsciiZHandler)
REGISTER_EXPORTER("hepevt_z", HepMC3HEPEVTZHandler)
REGISTER_EXPORTER("hepmc_lzma", HepMC3AsciiLZMAHandler)
REGISTER_EXPORTER("hepevt_lzma", HepMC3HEPEVTLZMAHandler)
REGISTER_EXPORTER("hepmc_bz2", HepMC3AsciiBZ2Handler)
REGISTER_EXPORTER("hepevt_bz2", HepMC3HEPEVTBZ2Handler)
#endif
#endif

#ifdef HEPMC3_ROOTIO
#include <HepMC3/WriterRoot.h>
#include <HepMC3/WriterRootTree.h>
typedef cepgen::HepMC3Handler<WriterRoot> HepMC3RootHandler;
typedef cepgen::HepMC3Handler<WriterRootTree> HepMC3RootTreeHandler;
REGISTER_EXPORTER("hepmc_root", HepMC3RootHandler)
REGISTER_EXPORTER("hepmc_root_tree", HepMC3RootTreeHandler)
#endif

#ifdef HEPMC3_EXTRA_PLUGINS
#include <ConvertExample/include/WriterDOT.h>
#include <ConvertExample/include/WriterRootTreeOPAL.h>
typedef cepgen::HepMC3Handler<WriterDOT> HepMC3DOTHandler;
typedef cepgen::HepMC3Handler<WriterRootTreeOPAL> HepMC3RootTreeOPALHandler;
REGISTER_EXPORTER("hepmc_dot", HepMC3DOTHandler)
REGISTER_EXPORTER("hepmc_root_tree_opal", HepMC3RootTreeOPALHandler)
#endif
