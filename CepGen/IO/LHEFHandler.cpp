#include "CepGen/IO/LHEFHandler.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"

#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

namespace cepgen
{
  namespace output
  {
    LHEFHandler::LHEFHandler( const char* filename ) :
      ExportHandler( ExportHandler::LHE )
#if defined ( HEPMC_LHEF )
      , lhe_output_( new LHEF::Writer( filename ) )
#elif defined ( PYTHIA_LHEF )
      , pythia_( new Pythia8::Pythia ), lhaevt_( new LHAevent )
#endif
    {
#if defined ( PYTHIA_LHEF )
      lhaevt_->openLHEF( filename );
#endif
    }

    LHEFHandler::~LHEFHandler()
    {
#if defined ( PYTHIA_LHEF )
      if ( lhaevt_ )
        lhaevt_->closeLHEF( false ); // we do not want to rewrite the init block
#endif
    }

    void
    LHEFHandler::initialise( const Parameters& params )
    {
      std::ostringstream oss_init;
      oss_init
        << "<!--\n"
        << "  ***** Sample generated with CepGen v" << version() << " *****\n"
        << "  * process: " << params.processName() << " (" << params.kinematics.mode << ")\n";
      if ( params.kinematics.mode != KinematicsMode::ElasticElastic ) {
        oss_init
          << "  * structure functions: " << params.kinematics.structure_functions->type << "\n";
        if ( !params.hadroniserName().empty() )
          oss_init
            << "  * hadroniser: " << params.hadroniserName() << "\n";
      }
      oss_init
        << "  *--- incoming state\n";
      if ( params.kinematics.cuts.initial.q2.valid() )
        oss_init
          << "  * Q2 range (GeV2): "
          << params.kinematics.cuts.initial.q2.min() << ", "
          << params.kinematics.cuts.initial.q2.max() << "\n";
      if ( params.kinematics.mode != KinematicsMode::ElasticElastic
        && params.kinematics.cuts.remnants.mass_single.valid() )
        oss_init
          << "  * remnants mass range (GeV/c2): "
          << params.kinematics.cuts.remnants.mass_single.min() << ", "
          << params.kinematics.cuts.remnants.mass_single.max() << "\n";
      oss_init
        << "  *--- central system\n";
      if ( params.kinematics.cuts.central.pt_single.valid() )
        oss_init
          << "  * single particle pt (GeV/c): "
          << params.kinematics.cuts.central.pt_single.min() << ", "
          << params.kinematics.cuts.central.pt_single.max() << "\n";
      if ( params.kinematics.cuts.central.energy_single.valid() )
        oss_init
          << "  * single particle energy (GeV): "
          << params.kinematics.cuts.central.energy_single.min() << ", "
          << params.kinematics.cuts.central.energy_single.max() << "\n";
      if ( params.kinematics.cuts.central.eta_single.valid() )
        oss_init
          << "  * single particle eta: "
          << params.kinematics.cuts.central.eta_single.min() << ", "
          << params.kinematics.cuts.central.eta_single.max() << "\n";
      if ( params.kinematics.cuts.central.pt_sum.valid() )
        oss_init
          << "  * total pt (GeV/c): "
          << params.kinematics.cuts.central.mass_sum.min() << ", "
          << params.kinematics.cuts.central.mass_sum.max() << "\n";
      if ( params.kinematics.cuts.central.mass_sum.valid() )
        oss_init
          << "  * total invariant mass (GeV/c2): "
          << params.kinematics.cuts.central.mass_sum.min() << ", "
          << params.kinematics.cuts.central.mass_sum.max() << "\n";
      oss_init
        << "  **************************************************\n"
        << "-->";
#if defined ( HEPMC_LHEF )
      lhe_output_->headerBlock() << oss_init.str();
      //params.dump( lhe_output_->initComments(), false );
      LHEF::HEPRUP run = lhe_output_->heprup;
      run.IDBMUP = { (int)params.kinematics.incoming_beams.first.pdg, (int)params.kinematics.incoming_beams.second.pdg };
      run.EBMUP = { (double)params.kinematics.incoming_beams.first.pz, (double)params.kinematics.incoming_beams.second.pz };
      run.NPRUP = 1;
      run.resize();
      run.XSECUP[0] = params.integrator.result;
      run.XERRUP[0] = params.integrator.err_result;
      run.XMAXUP[0] = 1.;
      run.LPRUP[0] = 1;
      lhe_output_->heprup = run;
      lhe_output_->init();
#elif defined ( PYTHIA_LHEF )
      oss_init << std::endl; // LHEF is usually not beautifully parsed as a standard XML...
      lhaevt_->addComments( oss_init.str() );
      lhaevt_->initialise( params );
      pythia_->settings.mode( "Beams:frameType", 5 );
      pythia_->settings.mode( "Next:numberCount", 0 ); // remove some of the Pythia output
      pythia_->settings.flag( "ProcessLevel:all", false ); // we do not want Pythia to interfere...
      pythia_->setLHAupPtr( lhaevt_.get() );
      pythia_->init();
      lhaevt_->initLHEF();
#endif
    }

    void
    LHEFHandler::operator<<( const Event& ev )
    {
#if defined ( HEPMC_LHEF )
      LHEF::HEPEUP out;
      out.heprup = &lhe_output_->heprup;
      out.XWGTUP = 1.;
      out.XPDWUP = std::pair<double,double>( 0., 0. );
      out.SCALUP = 0.;
      out.AQEDUP = constants::alphaEM;
      out.AQCDUP = constants::alphaQCD;
      out.NUP = ev.numParticles();
      out.resize();
      for ( unsigned short ip = 0; ip < ev.numParticles(); ++ip ) {
        const Particle part = ev.at( ip );
        out.IDUP[ip] = part.integerPdgId(); // PDG id
        out.ISTUP[ip] = (short)part.status(); // status code
        out.MOTHUP[ip] = std::pair<int,int>( ( part.mothers().size() > 0 ) ? *part.mothers().begin()+1 : 0, ( part.mothers().size() > 1 ) ? *part.mothers().rbegin()+1 : 0 ); // mothers
        out.ICOLUP[ip] = std::pair<int,int>( 0, 0 );
        out.PUP[ip] = std::vector<double>( { { part.momentum().px(), part.momentum().py(), part.momentum().pz(), part.energy(), part.mass() } } ); // momentum
        out.VTIMUP[ip] = 0.; // invariant lifetime
        out.SPINUP[ip] = 0.;
      }
      lhe_output_->eventComments() << "haha";
      lhe_output_->hepeup = out;
      lhe_output_->writeEvent();
#elif defined ( PYTHIA_LHEF )
      lhaevt_->feedEvent( 0, ev );
      pythia_->next();
      lhaevt_->eventLHEF();
#endif
    }

    void
    LHEFHandler::setCrossSection( double xsect, double xsect_err )
    {
#if defined ( PYTHIA_LHEF )
      lhaevt_->setCrossSection( 0, xsect, xsect_err );
#endif
    }

    //---------------------------------------------------------------------------------------------
    // Define LHA event record if one uses Pythia to store the LHE
    //---------------------------------------------------------------------------------------------

#if defined ( PYTHIA_LHEF )
    LHEFHandler::LHAevent::LHAevent() : LHAup( 3 )
    {}

    void
    LHEFHandler::LHAevent::initialise( const Parameters& params )
    {
      setBeamA( (short)params.kinematics.incoming_beams.first.pdg, params.kinematics.incoming_beams.first.pz );
      setBeamB( (short)params.kinematics.incoming_beams.second.pdg, params.kinematics.incoming_beams.second.pz );
      addProcess( 0, params.integrator.result, params.integrator.err_result, 100. );
    }

    void
    LHEFHandler::LHAevent::addComments( const std::string& comments )
    {
      osLHEF << comments;
    }

    void
    LHEFHandler::LHAevent::setCrossSection( unsigned short proc_id, double xsect, double xsect_err )
    {
      setXSec( proc_id, xsect );
      setXErr( proc_id, xsect_err );
    }

    void
    LHEFHandler::LHAevent::feedEvent( unsigned short proc_id, const Event& ev, bool full_event )
    {
      const double scale = ev.getOneByRole( Particle::Intermediate ).mass();
      setProcess( proc_id, 1., scale, constants::alphaEM, constants::alphaQCD );

      const Particle& ip1 = ev.getOneByRole( Particle::IncomingBeam1 ), &ip2 = ev.getOneByRole( Particle::IncomingBeam2 );
      const Particles& op1 = ev[Particle::OutgoingBeam1], &op2 = ev[Particle::OutgoingBeam2];
      const double q2_1 = -( ip1.momentum()-op1[0].momentum() ).mass2(), q2_2 = -( ip2.momentum()-op2[0].momentum() ).mass2();
      const double x1 = q2_1/( q2_1+op1[0].mass2()-ip1.mass2() ), x2 = q2_2/( q2_2+op2[0].mass2()-ip2.mass2() );
      setIdX( ip1.integerPdgId(), ip2.integerPdgId(), x1, x2 );

      short parton1_pdgid = 0, parton2_pdgid = 0;
      for ( const auto& part : ev.particles() ) {
        short pdg_id = part.integerPdgId(), status = 0, moth1 = 0, moth2 = 0;
        switch ( part.role() ) {
          case Particle::Parton1:
          case Particle::Parton2: {
            if ( part.role() == Particle::Parton1 )
              parton1_pdgid = part.integerPdgId();
            if ( part.role() == Particle::Parton2 )
              parton2_pdgid = part.integerPdgId();
            if ( !full_event )
              continue;
            status = -2; // conserving xbj/Q2
          } break;
          case Particle::Intermediate: {
            if ( !full_event )
              continue;
            status = 2;
            if ( pdg_id == 0 )
              pdg_id = ev.at( *part.mothers().begin() ).integerPdgId();
          } break;
          case Particle::IncomingBeam1:
          case Particle::IncomingBeam2: {
            if ( !full_event )
              continue;
            status = -9;
          } break;
          case Particle::OutgoingBeam1:
          case Particle::OutgoingBeam2:
          case Particle::CentralSystem: {
            status = (short)part.status();
            if ( status != 1 )
              continue;
          } break;
          default: break;
        }
        if ( full_event ) {
          const auto& mothers = part.mothers();
          if ( mothers.size() > 0 )
            moth1 = *mothers.begin()+1;
          if ( mothers.size() > 1 )
            moth2 = *mothers.rbegin()+1;
        }
        const Particle::Momentum& mom = part.momentum();
        addParticle( pdg_id, status, moth1, moth2, 0, 0, mom.px(), mom.py(), mom.pz(), mom.energy(), mom.mass(), 0. ,0., 0. );
      }
      setPdf( parton1_pdgid, parton2_pdgid, x1, x2, scale, 0., 0., true );
    }
#endif
  }
}
