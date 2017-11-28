#ifdef LIBHEPMC

#include "HepMCHandler.h"

#ifdef HEPMC_VERSION3
#include "LHEFHandler.h"

#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

namespace CepGen
{
  namespace OutputHandler
  {
    LHEFHandler::LHEFHandler( const char* filename ) :
      ExportHandler( ExportHandler::LHE ),
      lhe_output_( new LHEF::Writer( filename ) )
    {}

    void
    LHEFHandler::initialise( const Parameters& params )
    {
      lhe_output_->headerBlock()
        << "<!--\n"
        << "***** Sample generated with CepGen v" << version() << " *****\n"
        << "* process: " << params.processName() << " (" << params.kinematics.mode << ")\n"
        << "* structure functions: " << params.kinematics.structure_functions << "\n"
        << "*--- incoming state\n"
        << "* Q² range (GeV²): " << params.kinematics.cuts.initial.at( Cuts::q2 ) << "\n"
        << "* remnants mass range (GeV): " << params.kinematics.cuts.remnants.at( Cuts::mass ) << "\n"
        << "*--- central system\n"
        << "* single particle pT (GeV): " << params.kinematics.cuts.central.at( Cuts::pt_single ) << "\n"
        << "* single particle energy (GeV): " << params.kinematics.cuts.central.at( Cuts::energy_single ) << "\n"
        << "* single particle eta: " << params.kinematics.cuts.central.at( Cuts::eta_single ) << "\n"
        << "**************************************************\n"
        << "-->";
      //params.dump( lhe_output_->initComments(), false );
      LHEF::HEPRUP run = lhe_output_->heprup;
      run.IDBMUP = params.kinematics.inpdg;
      run.EBMUP = params.kinematics.inp;
      run.NPRUP = 1;
      run.resize();
      run.XSECUP[0] = cross_sect_;
      run.XERRUP[0] = cross_sect_err_;
      run.XMAXUP[0] = 1.;
      run.LPRUP[0] = 1;
      lhe_output_->heprup = run;
      lhe_output_->init();
    }

    void
    LHEFHandler::operator<<( const Event* ev )
    {
      LHEF::HEPEUP out;
      out.heprup = &lhe_output_->heprup;
      out.XWGTUP = 1.;
      out.XPDWUP = std::pair<double,double>( 0., 0. );
      out.SCALUP = 0.;
      out.AQEDUP = Constants::alphaEM;
      out.AQCDUP = Constants::alphaQCD;
      out.NUP = ev->numParticles();
      out.resize();
      for ( unsigned short ip=0; ip<ev->numParticles(); ip++ ) {
        const Particle part = ev->getConstById( ip );
        out.IDUP[ip] = part.integerPdgId(); // PDG id
        out.ISTUP[ip] = part.status(); // status code
        out.MOTHUP[ip] = std::pair<int,int>( ( part.mothers().size() > 0 ) ? *part.mothers().begin()+1 : 0, ( part.mothers().size() > 1 ) ? *part.mothers().rbegin()+1 : 0 ); // mothers
        out.ICOLUP[ip] = std::pair<int,int>( 0, 0 );
        out.PUP[ip] = std::vector<double>( { { part.momentum().px(), part.momentum().py(), part.momentum().pz(), part.energy(), part.mass() } } ); // momentum
        out.VTIMUP[ip] = 0.; // invariant lifetime
        out.SPINUP[ip] = 0.;
      }
      lhe_output_->eventComments() << "haha";
      lhe_output_->hepeup = out;
      lhe_output_->writeEvent();
    }
  }
}

#endif

#endif
