#include "CepGen/IO/LHEFHandler.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Parameters.h"
#include "CepGen/Version.h"
#include "CepGen/Event/Event.h"

namespace CepGen
{
  namespace OutputHandler
  {
    LHEFHandler::LHEFHandler( const char* filename ) :
      ExportHandler( ExportHandler::LHE ),
#ifdef HEPMC_LHEF
      lhe_output_( new LHEF::Writer( filename ) )
#else
      pythia_( new Pythia8::Pythia ),
      py_lhe_output_( new Pythia8::LHAupFromPYTHIA8( &pythia_->process, &pythia_->info ) )
#endif
    {
#ifndef HEPMC_LHEF
      py_lhe_output_->openLHEF( filename );
#endif
    }

    LHEFHandler::~LHEFHandler()
    {
#ifndef HEPMC_LHEF
      py_lhe_output_->closeLHEF( true );
#endif
    }

    void
    LHEFHandler::initialise( const Parameters& params )
    {
#ifdef HEPMC_LHEF
      lhe_output_->headerBlock()
        << "<!--\n"
        << "***** Sample generated with CepGen v" << version() << " *****\n"
        << "* process: " << params.processName() << " (" << params.kinematics.mode << ")\n"
        << "* structure functions: " << params.kinematics.structure_functions << "\n"
        << "*--- incoming state\n";
      if ( params.kinematics.cuts.initial.count( Cuts::q2 ) )
        lhe_output_->headerBlock()
          << "* Q² range (GeV²): " << params.kinematics.cuts.initial.at( Cuts::q2 ) << "\n";
      if ( params.kinematics.cuts.remnants.count( Cuts::mass ) )
        lhe_output_->headerBlock()
          << "* remnants mass range (GeV): " << params.kinematics.cuts.remnants.at( Cuts::mass ) << "\n";
      lhe_output_->headerBlock()
        << "*--- central system\n";
      if ( params.kinematics.cuts.central.count( Cuts::pt_single ) )
        lhe_output_->headerBlock()
          << "* single particle pT (GeV): " << params.kinematics.cuts.central.at( Cuts::pt_single ) << "\n";
      if ( params.kinematics.cuts.central.count( Cuts::energy_single ) )
        lhe_output_->headerBlock()
          << "* single particle energy (GeV): " << params.kinematics.cuts.central.at( Cuts::energy_single ) << "\n";
      if ( params.kinematics.cuts.central.count( Cuts::eta_single ) )
        lhe_output_->headerBlock()
          << "* single particle eta: " << params.kinematics.cuts.central.at( Cuts::eta_single ) << "\n";
      lhe_output_->headerBlock()
        << "**************************************************\n"
        << "-->";
      //params.dump( lhe_output_->initComments(), false );
      LHEF::HEPRUP run = lhe_output_->heprup;
      run.IDBMUP = params.kinematics.inpdg;
      run.EBMUP = params.kinematics.inp;
      run.NPRUP = 1;
      run.resize();
      run.XSECUP[0] = params.integrator.result;
      run.XERRUP[0] = params.integrator.err_result;
      run.XMAXUP[0] = 1.;
      run.LPRUP[0] = 1;
      lhe_output_->heprup = run;
      lhe_output_->init();
#else
      lhaevt_ = std::unique_ptr<LHAevent>( new LHAevent( params ) );
      lhaevt_->setCrossSection( params.integrator.result, params.integrator.err_result );
      pythia_->settings.mode( "Beams:frameType", 5 );
      pythia_->settings.flag( "ProcessLevel:all", false );
      pythia_->setLHAupPtr( lhaevt_.get() );
      pythia_->init();
      py_lhe_output_->setInit();
      py_lhe_output_->initLHEF();
#endif
    }

    void
    LHEFHandler::operator<<( const Event& ev )
    {
#ifdef HEPMC_LHEF
      LHEF::HEPEUP out;
      out.heprup = &lhe_output_->heprup;
      out.XWGTUP = 1.;
      out.XPDWUP = std::pair<double,double>( 0., 0. );
      out.SCALUP = 0.;
      out.AQEDUP = Constants::alphaEM;
      out.AQCDUP = Constants::alphaQCD;
      out.NUP = ev.numParticles();
      out.resize();
      for ( unsigned short ip = 0; ip < ev.numParticles(); ++ip ) {
        const Particle part = ev.getConstById( ip );
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
#else
      lhaevt_->feedEvent( ev );
      pythia_->next();
      py_lhe_output_->setEvent();
      py_lhe_output_->eventLHEF();
#endif
    }

    void
    LHEFHandler::setCrossSection( double xsect, double xsect_err )
    {
#ifndef HEPMC_LHEF
      lhaevt_->setCrossSection( xsect, xsect_err );
      py_lhe_output_->updateSigma();
#endif
    }

    //---------------------------------------------------------------------------------------------
    // Define LHA event record if one uses Pythia to store the LHE
    //---------------------------------------------------------------------------------------------

#ifndef HEPMC_LHEF
    LHEFHandler::LHAevent::LHAevent( const Parameters& params )
    {
      std::cout << params.integrator.result << "/////" << params.integrator.err_result <<std::endl;
      setBeamA( (short)params.kinematics.inpdg.first, params.kinematics.inp.first );
      setBeamB( (short)params.kinematics.inpdg.second, params.kinematics.inp.second );
      addProcess( 0, params.integrator.result, params.integrator.err_result, 100. );
    }

    void
    LHEFHandler::LHAevent::setCrossSection( double xsect, double xsect_err )
    {
      setXSec( 0, xsect );
      setXErr( 0, xsect_err );
    }

    void
    LHEFHandler::LHAevent::feedEvent( const Event& ev )
    {
      setProcess( 0, 1., ev.getOneByRole( Particle::Intermediate ).mass(), Constants::alphaEM, Constants::alphaQCD );
      for ( const auto& part : ev.particles() ) {
        short pdg_id = part.integerPdgId(), status = 0, moth1 = 0, moth2 = 0;
        switch ( part.role() ) {
          case Particle::Parton1: case Particle::Parton2:
            status = -1;
            break;
          case Particle::Intermediate:
            status = 2;
            if ( pdg_id == 0 )
              pdg_id = ev.getConstById( *part.mothers().begin() ).integerPdgId();
            break;
          case Particle::IncomingBeam1: case Particle::IncomingBeam2:
            status = -9;
            break;
          case Particle::OutgoingBeam1: case Particle::OutgoingBeam2:
          case Particle::CentralSystem:
            status = 1;
            break;
          default: break;
        }
        const auto& mothers = part.mothers();
        if ( mothers.size() > 0 )
          moth1 = *mothers.begin()+1;
        if ( mothers.size() > 1 )
          moth2 = *mothers.rbegin()+1;
        const Particle::Momentum& mom = part.momentum();
        addParticle( pdg_id, status, moth1, moth2, 0, 0, mom.px(), mom.py(), mom.pz(), mom.energy(), mom.mass(), 0. ,0., 0. );
      }
      //listEvent();
    }
#endif
  }
}

