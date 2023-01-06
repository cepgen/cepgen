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

#include <sstream>

#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Parameters.h"

using namespace std;  // account for improper scoping in following includes
#include <HepMC3/LHEF.h>

namespace cepgen {
  namespace io {
    /**
     * \brief Handler for the LHE file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class LHEFHepMCHandler : public ExportModule {
    public:
      /// Class constructor
      explicit LHEFHepMCHandler(const ParametersList&);

      static ParametersDescription description();

      void initialise() override;
      /// Writer operator
      void operator<<(const Event&) override;
      void setCrossSection(double, double) override;

    private:
      /// Writer object (from HepMC)
      std::unique_ptr<LHEF::Writer> lhe_output_;
      LHEF::HEPRUP run_;
      bool compress_;
    };

    LHEFHepMCHandler::LHEFHepMCHandler(const ParametersList& params)
        : ExportModule(params),
          lhe_output_(new LHEF::Writer(steer<std::string>("filename"))),
          compress_(steer<bool>("compress")) {}

    void LHEFHepMCHandler::setCrossSection(double cross_section, double err) {
      lhe_output_->heprup.NPRUP = 1;
      lhe_output_->heprup.resize();
      lhe_output_->heprup.XMAXUP[0] = 1.;
      lhe_output_->heprup.LPRUP[0] = 1;
      lhe_output_->heprup.XSECUP[0] = cross_section;
      lhe_output_->heprup.XERRUP[0] = err;
    }

    void LHEFHepMCHandler::initialise() {
      lhe_output_->headerBlock() << "<!--\n" << banner() << "\n-->";
      //--- first specify information about the run
      lhe_output_->heprup.IDBMUP = {(int)rt_params_->kinematics().incomingBeams().positive().pdgId(),
                                    (int)rt_params_->kinematics().incomingBeams().negative().pdgId()};
      lhe_output_->heprup.EBMUP = {(double)rt_params_->kinematics().incomingBeams().positive().momentum().pz(),
                                   (double)rt_params_->kinematics().incomingBeams().negative().momentum().pz()};
      //--- ensure everything is properly parsed
      lhe_output_->init();
    }

    void LHEFHepMCHandler::operator<<(const Event& ev) {
      LHEF::HEPEUP out;
      out.heprup = &lhe_output_->heprup;
      out.XWGTUP = 1.;
      out.XPDWUP = std::pair<double, double>(0., 0.);
      out.SCALUP = 0.;
      out.AQEDUP = constants::ALPHA_EM;
      out.AQCDUP = constants::ALPHA_QCD;
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
      //lhe_output_->eventComments() << "haha";
      lhe_output_->hepeup = out;
      lhe_output_->writeEvent();
    }

    ParametersDescription LHEFHepMCHandler::description() {
      auto desc = ExportModule::description();
      desc.setDescription("HepMC 3-based LHEF output module");
      desc.add<std::string>("filename", "output.lhe").setDescription("Output filename");
      desc.add<bool>("compress", true);
      return desc;
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("lhef_hepmc", LHEFHepMCHandler)
