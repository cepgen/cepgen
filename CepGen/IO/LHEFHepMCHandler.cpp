#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"

#include <sstream>

#if !defined( HEPMC3 )
#  include "HepMC/Version.h"
#  ifndef HEPMC_VERSION_CODE // HepMC v2
#    error "HepMC v3 is required for the LHEF export!"
#  else // HepMC v3+
#    include "HepMC/LHEF.h"
#  endif
#else // HepMC v3+
#  include "HepMC3/LHEF.h"
#endif

namespace cepgen
{
  namespace output
  {
    /**
     * \brief Handler for the LHE file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class LHEFHepMCHandler : public GenericExportHandler
    {
      public:
        /// Class constructor
        explicit LHEFHepMCHandler( const ParametersList& );
        ~LHEFHepMCHandler() override;

        void initialise( const Parameters& ) override;
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override {}

      private:
        /// Writer object (from HepMC)
        std::unique_ptr<LHEF::Writer> lhe_output_;
        LHEF::HEPRUP run_;
    };

    LHEFHepMCHandler::LHEFHepMCHandler( const ParametersList& params ) :
      GenericExportHandler( GenericExportHandler::LHE ),
      lhe_output_( new LHEF::Writer( params.get<std::string>( "filename", "output.lhe" ) ) )
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
      out.NUP = ev.numParticles();
      out.resize();
      for ( unsigned short ip = 0; ip < ev.numParticles(); ++ip ) {
        const Particle part = ev[ip];
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

REGISTER_IO_MODULE( lhef, LHEFHepMCHandler )
