#include "CepGen/Modules/ExportModule.h"
#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Parameters.h"

#include <sstream>

#ifdef HEPMC3
using namespace std; // account for improper scoping in following includes
#  include "HepMC3/LHEF.h"
#else
#  include "HepMC/Version.h"
#  ifdef HEPMC_VERSION_CODE // HepMC v3+
#    include "HepMC/LHEF.h"
#  else
#    define NO_LHEF
#  endif
#endif
#ifndef NO_LHEF

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Handler for the LHE file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class LHEFHepMCHandler : public ExportModule
    {
      public:
        /// Class constructor
        explicit LHEFHepMCHandler( const ParametersList& );

        void initialise( const Parameters& ) override;
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override {}

      private:
        /// Writer object (from HepMC)
        std::unique_ptr<LHEF::Writer> lhe_output_;
        LHEF::HEPRUP run_;
        bool compress_;
    };

    LHEFHepMCHandler::LHEFHepMCHandler( const ParametersList& params ) :
      ExportModule( params ),
      lhe_output_( new LHEF::Writer( params.get<std::string>( "filename", "output.lhe" ) ) ),
      compress_( params.get<bool>( "compress", true ) )
    {}

    void
    LHEFHepMCHandler::initialise( const Parameters& params )
    {
      lhe_output_->headerBlock()
        << "<!--\n" << banner( params ) << "\n-->";
      //--- first specify information about the run
      LHEF::HEPRUP run = lhe_output_->heprup;
      run.IDBMUP = { (int)params.kinematics.incoming_beams.first.pdg, (int)params.kinematics.incoming_beams.second.pdg };
      run.EBMUP = { (double)params.kinematics.incoming_beams.first.pz, (double)params.kinematics.incoming_beams.second.pz };
      run.NPRUP = 1;
      run.resize();
      run.XSECUP[0] = params.integration().result;
      run.XERRUP[0] = params.integration().err_result;
      run.XMAXUP[0] = 1.;
      run.LPRUP[0] = 1;
      lhe_output_->heprup = run;
      //--- ensure everything is properly parsed
      lhe_output_->init();
    }

    void
    LHEFHepMCHandler::operator<<( const Event& ev )
    {
      LHEF::HEPEUP out;
      out.heprup = &lhe_output_->heprup;
      out.XWGTUP = 1.;
      out.XPDWUP = std::pair<double,double>( 0., 0. );
      out.SCALUP = 0.;
      out.AQEDUP = constants::ALPHA_EM;
      out.AQCDUP = constants::ALPHA_QCD;
      const auto& particles = compress_
        ? ev.compress().particles()
        : ev.particles();
      out.NUP = particles.size();
      out.resize();
      for ( unsigned short ip = 0; ip < particles.size(); ++ip ) {
        const Particle& part = particles[ip];
        out.IDUP[ip] = part.integerPdgId(); // PDG id
        out.ISTUP[ip] = (short)part.status(); // status code
        out.PUP[ip] = part.momentum().pVector(); // momentum
        out.MOTHUP[ip] = { // mothers
          part.mothers().size() > 0 ? *part.mothers(). begin()+1 : 0,
          part.mothers().size() > 1 ? *part.mothers().rbegin()+1 : 0
        };
        out.ICOLUP[ip] = { 0, 0 };
        out.VTIMUP[ip] = 0.; // invariant lifetime
        out.SPINUP[ip] = 0.;
      }
      //lhe_output_->eventComments() << "haha";
      lhe_output_->hepeup = out;
      lhe_output_->writeEvent();
    }
  }
}

REGISTER_IO_MODULE( "lhef", LHEFHepMCHandler )
#endif
