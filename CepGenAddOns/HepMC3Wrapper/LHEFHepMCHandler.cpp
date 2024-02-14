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
      desc.add<bool>("compress", true);
      return desc;
    }

    bool operator<<(const Event& ev) override {
      LHEF::HEPEUP out;
      out.heprup = &lhe_output_->heprup;
      out.XWGTUP = 1.;
      out.XPDWUP = std::pair<double, double>(0., 0.);
      out.SCALUP = 0.;
      out.AQEDUP = ev.metadata("alphaEM");
      out.AQCDUP = ev.metadata("alphaS");
      const auto& particles = compress_ ? ev.compress().particles() : ev.particles();
      out.NUP = particles.size();
      out.resize();
      for (unsigned short ip = 0; ip < particles.size(); ++ip) {
        const Particle& part = particles[ip];
        out.IDUP[ip] = part.integerPdgId();    // PDG id
        out.ISTUP[ip] = (short)part.status();  // status code
        std::copy(part.momentum().pVector().begin(), part.momentum().pVector().end(),
                  out.PUP[ip].begin());  // momentum
        out.MOTHUP[ip] = {               // mothers
                          part.mothers().size() > 0 ? *part.mothers().begin() + 1 : 0,
                          part.mothers().size() > 1 ? *part.mothers().rbegin() + 1 : 0};
        out.ICOLUP[ip] = {0, 0};
        out.VTIMUP[ip] = 0.;  // invariant lifetime
        out.SPINUP[ip] = 0.;
      }
      lhe_output_->hepeup = out;
      lhe_output_->writeEvent();
    }

    void setCrossSection(const Value& cross_section) override {
      lhe_output_->heprup.NPRUP = 1;
      lhe_output_->heprup.resize();
      lhe_output_->heprup.XMAXUP[0] = 1.;
      lhe_output_->heprup.LPRUP[0] = 1;
      lhe_output_->heprup.XSECUP[0] = (double)cross_section;
      lhe_output_->heprup.XERRUP[0] = cross_section.uncertainty();
    }

  private:
    void initialise() override {
      lhe_output_->headerBlock() << "<!--\n" << banner() << "\n-->";
      // run information
      lhe_output_->heprup.IDBMUP = {(int)runParameters().kinematics().incomingBeams().positive().pdgId(),
                                    (int)runParameters().kinematics().incomingBeams().negative().pdgId()};
      lhe_output_->heprup.EBMUP = {(double)runParameters().kinematics().incomingBeams().positive().momentum().pz(),
                                   (double)runParameters().kinematics().incomingBeams().negative().momentum().pz()};
      lhe_output_->init();  // ensure everything is properly parsed
    }

    const std::unique_ptr<LHEF::Writer> lhe_output_;
    const bool compress_;
  };
}  // namespace cepgen
REGISTER_EXPORTER("lhef_hepmc", LHEFHepMCHandler);
