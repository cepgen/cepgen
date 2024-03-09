/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2024  Laurent Forthomme
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

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Utils/Value.h"

using namespace std;  // account for improper scoping in following include
#include <HepMC3/LHEF.h>

namespace cepgen {
  /// Handler for the LHE file output
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Sep 2016
  class LHEFHepMCHandler : public EventExporter {
  public:
    explicit LHEFHepMCHandler(const ParametersList& params)
        : EventExporter(params),
          lhe_output_(new LHEF::Writer(steer<std::string>("filename"))),
          compress_(steer<bool>("compress")) {}

    static ParametersDescription description() {
      auto desc = EventExporter::description();
      desc.setDescription("HepMC 3-based LHEF output module");
      desc.add<std::string>("filename", "output.lhe").setDescription("Output filename");
      desc.add<bool>("compress", false);
      return desc;
    }

    bool operator<<(const Event& cg_ev) override {
      if (!header_initialised_) {
        lhe_output_->init();  // ensure everything is properly parsed
        header_initialised_ = true;
      }
      auto& hepeup = lhe_output_->hepeup;
      hepeup.heprup = &lhe_output_->heprup;
      hepeup.XWGTUP = 1.;
      hepeup.XPDWUP = std::pair<double, double>(0., 0.);
      hepeup.SCALUP = 0.;
      hepeup.AQEDUP = cg_ev.metadata("alphaEM");
      hepeup.AQCDUP = cg_ev.metadata("alphaS");
      const auto cg_particles = compress_ ? cg_ev.compress().particles() : cg_ev.particles();
      hepeup.resize(cg_particles.size());
      for (unsigned short ip = 0; ip < hepeup.NUP; ++ip) {
        const auto& cg_part = cg_particles.at(ip);
        hepeup.IDUP[ip] = cg_part.integerPdgId();    // PDG id
        hepeup.ISTUP[ip] = (short)cg_part.status();  // status code
        hepeup.PUP[ip] = {cg_part.momentum().px(),
                          cg_part.momentum().py(),
                          cg_part.momentum().pz(),
                          cg_part.momentum().energy(),
                          cg_part.momentum().mass()};
        hepeup.MOTHUP[ip] = {// mothers
                             cg_part.mothers().size() > 0 ? *cg_part.mothers().begin() + 1 : 0,
                             cg_part.mothers().size() > 1 ? *cg_part.mothers().rbegin() + 1 : 0};
        hepeup.ICOLUP[ip] = {0, 0};
        hepeup.VTIMUP[ip] = 0.;  // invariant lifetime
        hepeup.SPINUP[ip] = 0.;
      }
      lhe_output_->writeEvent();
      return true;
    }

    void setCrossSection(const Value& cross_section) override {
      lhe_output_->heprup.XSECUP[0] = (double)cross_section;
      lhe_output_->heprup.XERRUP[0] = cross_section.uncertainty();
    }

  private:
    void initialise() override {
      lhe_output_->headerBlock() << "<!--\n" << banner() << "\n-->";
      if (runParameters().hasProcess()) {  // run information only specified if process (and kinematics) is specified
        lhe_output_->heprup.IDBMUP = {(int)runParameters().kinematics().incomingBeams().positive().integerPdgId(),
                                      (int)runParameters().kinematics().incomingBeams().negative().integerPdgId()};
        lhe_output_->heprup.EBMUP = {(double)runParameters().kinematics().incomingBeams().positive().momentum().pz(),
                                     (double)runParameters().kinematics().incomingBeams().negative().momentum().pz()};
      }
      lhe_output_->heprup.resize(1);
      lhe_output_->heprup.XMAXUP[0] = 1.;
      lhe_output_->heprup.LPRUP[0] = 1;
      lhe_output_->heprup.XSECUP[0] = 0.;  // placeholders
      lhe_output_->heprup.XERRUP[0] = 0.;
    }

    const std::unique_ptr<LHEF::Writer> lhe_output_;
    const bool compress_;
    bool header_initialised_{false};
  };
}  // namespace cepgen
REGISTER_EXPORTER("lhef_hepmc", LHEFHepMCHandler);
