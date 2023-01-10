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

#include <HepMC/Version.h>

#include <memory>

#include "CepGen/Core/EventExporter.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Parameters.h"
#include "CepGenAddOns/HepMC2Wrapper/HepMC2EventInterface.h"

using namespace HepMC;

namespace cepgen {
  /// Handler for the HepMC file output
  /// \tparam T HepMC writer handler (format-dependent)
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Sep 2016
  template <typename T>
  class HepMC2Handler : public EventExporter {
  public:
    /// Class constructor
    explicit HepMC2Handler(const ParametersList&);
    ~HepMC2Handler();

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
  };

  template <typename T>
  HepMC2Handler<T>::HepMC2Handler(const ParametersList& params)
      : EventExporter(params), output_(new T(steer<std::string>("filename").c_str())), xs_(new GenCrossSection) {
    CG_INFO("HepMC") << "Interfacing module initialised "
                     << "for HepMC version " << HEPMC_VERSION << ".";
  }

  template <typename T>
  HepMC2Handler<T>::~HepMC2Handler() {}

  template <typename T>
  void HepMC2Handler<T>::operator<<(const Event& evt) {
    CepGenEvent event(evt);
    // general information
    event.set_cross_section(*xs_);
    event.set_event_number(event_num_++);
    output_->write_event(&event);
  }

  template <typename T>
  void HepMC2Handler<T>::setCrossSection(double cross_section, double cross_section_err) {
    xs_->set_cross_section(cross_section, cross_section_err);
  }

  template <typename T>
  ParametersDescription HepMC2Handler<T>::description() {
    auto desc = EventExporter::description();
    desc.setDescription("HepMC2 ASCII file output module");
    desc.add<std::string>("filename", "output.hepmc").setDescription("Output filename");
    return desc;
  }
}  // namespace cepgen

//----------------------------------------------------------------------
// Defining the various templated plugins made available by this
// specific version of HepMC
//----------------------------------------------------------------------

//--- HepMC version 2 and below
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/IO_GenEvent.h"
typedef cepgen::HepMC2Handler<IO_GenEvent> HepMC2GenEventHandler;
typedef cepgen::HepMC2Handler<IO_AsciiParticles> HepMC2AsciiHandler;
REGISTER_EXPORTER("hepmc2", HepMC2GenEventHandler)
REGISTER_EXPORTER("hepmc2_ascii", HepMC2AsciiHandler)
