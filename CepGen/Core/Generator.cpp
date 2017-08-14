#include "CepGen/Generator.h"

namespace CepGen
{
  Generator::Generator() :
    cross_section_( -1. ), cross_section_error_( -1. ), has_cross_section_( false )
  {
    Debugging( "Generator initialized" );
    try { printHeader(); } catch ( Exception& e ) { e.dump(); }
    srand( time( 0 ) ); // Random number initialization
    this->parameters = std::unique_ptr<Parameters>( new Parameters );
  }

  Generator::Generator( Parameters* ip ) :
    parameters( ip )
  {}

  Generator::~Generator()
  {
    if ( parameters->generation && parameters->process() && parameters->process()->numGeneratedEvents()>0 ) {
      Information( Form( "Mean generation time / event: %.3f ms", parameters->process()->totalGenerationTime()*1.e3/parameters->process()->numGeneratedEvents() ) );
    }
  }

  void
  Generator::clearRun()
  {
    parameters->vegas.first_run = true;
    has_cross_section_ = false; // force the recreation of the Vegas instance
    cross_section_ = cross_section_error_ = -1.;
  }

  void
  Generator::setParameters( Parameters& ip )
  {
    parameters = std::unique_ptr<Parameters>( new Parameters( ip ) ); // copy constructor
  }

  void
  Generator::printHeader()
  {
    std::string tmp;
    std::ostringstream os; os << std::endl;
    std::ifstream hf( "README" );
    if ( !hf.good() ) throw Exception( __PRETTY_FUNCTION__, "Failed to open README file", JustWarning );
    while ( true ) {
      if ( !hf.good() ) break;
      getline( hf, tmp );
      os << "\n " << tmp;
    }
    hf.close();
    Information( os.str().c_str() );
  }

  void
  Generator::computeXsection( double& xsec, double& err )
  {
    // first destroy and recreate the Vegas instance
    if ( !vegas_ ) {
      vegas_ = std::unique_ptr<Vegas>( new Vegas( numDimensions(), f, parameters.get() ) );
    }
    else if ( vegas_->dimensions() != numDimensions() ) {
      vegas_.reset( new Vegas( numDimensions(), f, parameters.get() ) );
    }

    if ( Logger::get().level>=Logger::Debug ) {
      std::ostringstream topo; topo << parameters->process_mode;
      Debugging( Form( "New Vegas instance created\n\t"
                       "Considered topology: %s case\n\t"
                       "Will proceed with %d-dimensional integration", topo.str().c_str(), numDimensions() ) );
    }

    Information( "Starting the computation of the process cross-section" );

    try { prepareFunction(); } catch ( Exception& e ) { e.dump(); }

    has_cross_section_ = ( vegas_->integrate( cross_section_, cross_section_error_ ) == 0 );

    xsec = cross_section_;
    err = cross_section_error_;

    Information( Form( "Total cross section: %f +/- %f pb", xsec, err ) );
  }

  Event*
  Generator::generateOneEvent()
  {
    bool good = false;
    if ( !has_cross_section_ ) {
      computeXsection( cross_section_, cross_section_error_ );
    }
    while ( !good ) { good = vegas_->generateOneEvent(); }

    last_event = this->parameters->last_event;
    return last_event.get();
  }

  void
  Generator::prepareFunction()
  {
    if ( !parameters->process() ) {
      throw Exception( __PRETTY_FUNCTION__, "No process defined!", FatalError );
    }
    Kinematics kin = parameters->kinematics;
    parameters->process()->addEventContent();
    parameters->process()->setKinematics( kin );
    Debugging( "Function prepared to be integrated!" );
  }

  double
  f( double* x, size_t ndim, void* params )
  {
    Timer tmr;
    bool hadronised;
    double num_hadr_trials;
    std::ostringstream os;

    Parameters* p = static_cast<Parameters*>( params );

    //float now = tmr.elapsed();
    const Particle::Momentum p1( 0., 0.,  p->kinematics.in1p ),
                             p2( 0., 0., -p->kinematics.in2p );
    p->process()->setIncomingKinematics( p1, p2 ); // at some point introduce non head-on colliding beams?
    //PrintMessage( Form( "0 - after setting the kinematics: %.3e", tmr.elapsed()-now ) ); now = tmr.elapsed();

    DebuggingInsideLoop( Form( "Function f called -- some parameters:\n\t"
                               "  pz(p1) = %5.2f  pz(p2) = %5.2f\n\t"
                               "  remnant mode: %d",
                               p->kinematics.in1p, p->kinematics.in2p, p->remnant_mode ) );

    p->process()->clearEvent();
    //PrintMessage( Form( "1 - after clearing the event: %.3e", tmr.elapsed()-now ) ); now = tmr.elapsed();
  
    Event* ev = p->process()->event();

    //PrintMessage( Form( "2: %.3e", tmr.elapsed()-now ) ); now = tmr.elapsed();

    if ( p->vegas.first_run ) {
      //PrintMessage( Form(  "3: %.3e", tmr.elapsed()-now ) ); now = tmr.elapsed();

      if ( Logger::get().level >= Logger::Debug ) {
        std::ostringstream oss; oss << p->process_mode;
        Debugging( Form( "Process mode considered: %s", oss.str().c_str() ) );
      }

      //--- add outgoing protons or remnants
      switch ( p->process_mode ) {
        case Kinematics::ElectronProton: default: { InError( "Not handled yet!" ); }
        case Kinematics::ElasticElastic: break; // nothing to change in the event
        case Kinematics::ElasticInelastic:
        case Kinematics::InelasticElastic: // set one of the outgoing protons to be fragmented
          ev->getOneByRole( Particle::OutgoingBeam1 )->setPdgId( Particle::uQuark ); break;
        case Kinematics::InelasticInelastic: // set both the outgoing protons to be fragmented
          ev->getOneByRole( Particle::OutgoingBeam1 )->setPdgId( Particle::uQuark );
          ev->getOneByRole( Particle::OutgoingBeam2 )->setPdgId( Particle::uQuark );
          break;
      }
      //PrintMessage( Form( "4 - after preparing the event kinematics: %.3e", tmr.elapsed()-now ) ); now = tmr.elapsed();

      //--- prepare the function to be integrated
      p->process()->prepareKinematics();

      //--- add outgoing leptons
      Particle* out1 = ev->getOneByRole( Particle::CentralParticle1 ),
               *out2 = ev->getOneByRole( Particle::CentralParticle2 );
      out1->setPdgId( p->kinematics.pair ); out1->setMass( Particle::massFromPDGId( p->kinematics.pair ) );
      out2->setPdgId( p->kinematics.pair ); out2->setMass( Particle::massFromPDGId( p->kinematics.pair ) );

      p->process()->clearRun();
      p->vegas.first_run = false;
    }

    p->process()->setPoint( ndim, x );
    if ( Logger::get().level >= Logger::DebugInsideLoop ) {
      std::ostringstream oss; for ( unsigned int i=0; i<ndim; i++ ) { oss << x[i] << " "; }
      DebuggingInsideLoop( Form( "Computing dim-%d point ( %s)", ndim, oss.str().c_str() ) );
    }

    //--- from this step on, the phase space point is supposed to be set by 'GenericProcess::setPoint()'

    p->process()->beforeComputeWeight();

    tmr.reset();
    //PrintMessage( Form( "5: %.3e", tmr.elapsed()-now ) ); now = tmr.elapsed();
    double ff = p->process()->computeWeight();
    //PrintMessage( Form( "6: %.3e", tmr.elapsed()-now ) ); now = tmr.elapsed();

    if ( ff<0. ) return 0.;

    if ( p->store ) { // MC events generation
      p->process()->fillKinematics( false );

      ev->time_generation = tmr.elapsed();

      if ( p->hadroniser() && p->process_mode!=Kinematics::ElasticElastic ) {

        Debugging( Form( "Event before calling the hadroniser (%s)", p->hadroniser()->name().c_str() ) );
        if ( Logger::get().level>=Logger::Debug ) ev->dump();

        num_hadr_trials = 0;
        do {
          try { hadronised = p->hadroniser()->hadronise( ev ); } catch ( Exception& e ) { e.dump(); }

          if ( num_hadr_trials>0 ) { Debugging( Form( "Hadronisation failed. Trying for the %dth time", num_hadr_trials+1 ) ); }

          num_hadr_trials++;
        } while ( !hadronised && num_hadr_trials<=p->hadroniser_max_trials );
        if ( !hadronised ) return 0.; //FIXME

        ev->num_hadronisation_trials = num_hadr_trials;

        Debugging( Form( "Event hadronisation succeeded after %d trial(s)", ev->num_hadronisation_trials ) );

        if ( num_hadr_trials>p->hadroniser_max_trials ) return 0.; //FIXME

        Debugging( Form( "Event after calling the hadroniser (%s)", p->hadroniser()->name().c_str() ) );
        if ( Logger::get().level>=Logger::Debug ) ev->dump();
      }
      ev->time_total = tmr.elapsed();
      p->process()->addGenerationTime( ev->time_total );

      Debugging( Form( "Generation time:       %5.6f sec\n\t"
                       "Total time (gen+hadr): %5.6f sec",
                       ev->time_generation,
                       ev->time_total ) );

      *(p->last_event) = *( ev );
      //ev->Store(p->file);
    }

    if ( Logger::get().level>=Logger::DebugInsideLoop ) {
      os.str( "" ); for ( unsigned int i=0; i<ndim; i++ ) { os << Form( "%10.8f ", x[i] ); }
      Debugging( Form( "f value for dim-%d point ( %s): %4.4e", ndim, os.str().c_str(), ff ) );
    }

    return ff;
  }
}

