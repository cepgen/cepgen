#include "CepGen/Core/ExportModule.h"
#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Parameters.h"

#include <sstream>

#ifdef HEPMC3
using namespace std;  // account for improper scoping in following includes
#include "HepMC3/LHEF.h"
#else
#include "HepMC/Version.h"
#ifdef HEPMC_VERSION_CODE  // HepMC v3+
#include "HepMC/LHEF.h"
#else
#define NO_LHEF
#endif
#endif
#ifndef NO_LHEF

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
      static std::string description() { return "HepMC 3-based LHEF output module"; }

      void initialise(const Parameters&) override;
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
          lhe_output_(new LHEF::Writer(params.get<std::string>("filename", "output.lhe"))),
          compress_(params.get<bool>("compress", true)) {}

    void LHEFHepMCHandler::setCrossSection(double cross_section, double err) {
      lhe_output_->heprup.NPRUP = 1;
      lhe_output_->heprup.resize();
      lhe_output_->heprup.XMAXUP[0] = 1.;
      lhe_output_->heprup.LPRUP[0] = 1;
      lhe_output_->heprup.XSECUP[0] = cross_section;
      lhe_output_->heprup.XERRUP[0] = err;
    }

    void LHEFHepMCHandler::initialise(const Parameters& params) {
      lhe_output_->headerBlock() << "<!--\n" << banner(params) << "\n-->";
      //--- first specify information about the run
      lhe_output_->heprup.IDBMUP = {(int)params.kinematics.incoming_beams.first.pdg,
                                    (int)params.kinematics.incoming_beams.second.pdg};
      lhe_output_->heprup.EBMUP = {(double)params.kinematics.incoming_beams.first.pz,
                                   (double)params.kinematics.incoming_beams.second.pz};
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
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("lhef_hepmc", LHEFHepMCHandler)
#endif
