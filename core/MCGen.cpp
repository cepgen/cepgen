#include "MCGen.h"

MCGen::MCGen() :
  vegas_( 0 ), cross_section_( -1. ), cross_section_error_( -1. ), has_cross_section_( false )
{
  Debugging( "Generator initialized" );
  
  try { printHeader(); } catch ( Exception& e ) { e.dump(); }
  
  srand( time( 0 ) ); // Random number initialization
  
  this->parameters = new Parameters;
}

MCGen::MCGen( Parameters *ip_ ) :
  parameters( ip_ ), vegas_( 0 )
{}

MCGen::~MCGen()
{
  if ( parameters->generation and parameters->process and parameters->process->numGeneratedEvents()>0 ) {
    Information( Form( "Mean generation time / event: %.3f ms", parameters->process->totalGenerationTime()*1.e3/parameters->process->numGeneratedEvents() ) );
  }

  if ( vegas_ ) delete vegas_;
  if ( parameters ) delete parameters;
}

void
MCGen::printHeader()
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
MCGen::buildVegas()
{
  if ( Logger::GetInstance()->Level>=Logger::Debug ) {
    std::ostringstream topo; topo << parameters->process_mode;
    Debugging( Form( "Considered topology: %s case\n\t"
		     "Will proceed with %d-dimensional integration", topo.str().c_str(), numDimensions() ) );
  }
  
  if ( vegas_ ) delete vegas_;
  vegas_ = new Vegas( numDimensions(), f, parameters );
}

void
MCGen::computeXsection( double& xsec, double& err )
{
  if ( !vegas_ ) buildVegas();

  Information( "Starting the computation of the process cross-section" );

  try { prepareFunction(); } catch ( Exception& e ) { e.dump(); }
  vegas_->integrate( xsec, err );
  
  cross_section_ = xsec;
  cross_section_error_ = err;
  has_cross_section_ = true;
  
  Information( Form( "Total cross section: %f +/- %f pb", xsec, err ) );
}

Event*
MCGen::generateOneEvent()
{
  bool good = false;
  if ( !has_cross_section_ ) {
    double xsec, err;
    computeXsection( xsec, err );
  }
  while ( !good ) { good = vegas_->generateOneEvent(); }

  last_event = this->parameters->last_event;
  return static_cast<Event*>( last_event );
}

void
MCGen::prepareFunction()
{
  if ( !parameters->process ) {
    throw Exception( __PRETTY_FUNCTION__, "No process defined!", FatalError );
  }
  Kinematics kin;
  kin.kinematics = static_cast<Kinematics::ProcessMode>( parameters->process_mode );
  /*kin.q1tmin = kin.q2tmin = 0.;
  kin.q1tmax = kin.q2tmax = 50.;*/
  kin.q2min = parameters->minq2;
  kin.q2max = parameters->maxq2;
  kin.qtmin = parameters->minqt;
  kin.qtmax = parameters->maxqt;
  kin.mode = parameters->mcut;
  kin.ptmin = parameters->minpt;
  kin.ptmax = parameters->maxpt;
  kin.ptdiffmin = parameters->minptdiff;
  kin.ptdiffmax = parameters->maxptdiff;
  kin.etamin = parameters->mineta;
  kin.etamax = parameters->maxeta;
  kin.massmin = parameters->minmass;
  kin.massmax = parameters->maxmass;
  kin.emin = parameters->minenergy;
  kin.emax = parameters->maxenergy;
  kin.mxmin = parameters->minmx;
  kin.mxmax = parameters->maxmx;
  kin.remnant_mode = parameters->remnant_mode;
  parameters->process->addEventContent();
  parameters->process->setKinematics( kin );
  Debugging( "Function prepared to be integrated!" );
}

double f( double* x, size_t ndim, void* params )
{
  Timer tmr;
  bool hadronised;
  double num_hadr_trials;
  std::ostringstream os;

  Parameters* p = static_cast<Parameters*>( params );

  //FIXME at some point introduce non head-on colliding beams ?

  const Particle::Momentum p1( 0., 0.,  p->in1p ),
                           p2( 0., 0., -p->in2p );
  p->process->setIncomingKinematics( p1, p2 );
  p->process->setPoint( ndim, x );

  if ( Logger::GetInstance()->Level>=Logger::DebugInsideLoop ) {
    os.str(""); for ( unsigned int i=0; i<ndim; i++ ) { os << x[i] << " "; }
    DebuggingInsideLoop( Form( "Computing dim-%d point ( %s)", ndim, os.str().c_str() ) );
  }

  tmr.reset();

  double ff = 0.;

  DebuggingInsideLoop( Form( "Function f called -- some parameters:\n\t"
                             "  pz(p1) = %5.2f  pz(p2) = %5.2f\n\t"
                             "  remnant mode: %d",
                             p->in1p, p->in2p, p->remnant_mode ) );
    
  //float now = tmr.elapsed();
  //std::cout << "0: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
  p->process->clearEvent();
  //std::cout << "1: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
  
  Event* ev = p->process->event();
  
  //std::cout << "2: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();*/

  if ( p->first_run ) {
    //std::cout << "3: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
    // Then add outgoing protons or remnants
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
    
    //std::cout << "4: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
    // Prepare the function to be integrated
    p->process->prepareKinematics();

    // Then add outgoing leptons
    Particle* out1 = ev->getOneByRole( Particle::CentralParticle1 ),
             *out2 = ev->getOneByRole( Particle::CentralParticle2 );
    out1->setPdgId( p->pair ); out1->setMass( Particle::massFromPDGId( p->pair ) );
    out2->setPdgId( p->pair ); out2->setMass( Particle::massFromPDGId( p->pair ) );

    p->process->clearRun();
    p->first_run = false;
  }
   
  p->process->beforeComputeWeight();

  //std::cout << "5: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
  ff = p->process->computeWeight();
  //std::cout << "6: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
  if (ff<0.) return 0.;
  
  if (p->store) { // MC events generation
    p->process->fillKinematics( false );
    
    ev->time_generation = tmr.elapsed();

    if ( p->hadroniser and p->process_mode!=Kinematics::ElasticElastic ) {
      
      Debugging( Form( "Event before calling the hadroniser (%s)", p->hadroniser->name().c_str() ) );
      if ( Logger::GetInstance()->Level>=Logger::Debug ) ev->dump();
      
      num_hadr_trials = 0;
      do {
        try { hadronised = p->hadroniser->hadronise( ev ); } catch ( Exception& e ) { e.dump(); }

        if ( num_hadr_trials>0 ) { Debugging( Form( "Hadronisation failed. Trying for the %dth time", num_hadr_trials+1 ) ); }
        
        num_hadr_trials++;
      } while ( !hadronised and num_hadr_trials<=p->hadroniser_max_trials );
      if ( !hadronised ) return 0.; //FIXME
      
      ev->num_hadronisation_trials = num_hadr_trials;

      Debugging( Form( "Event hadronisation succeeded after %d trial(s)", ev->num_hadronisation_trials ) );

      if ( num_hadr_trials>p->hadroniser_max_trials ) return 0.; //FIXME
      
      Debugging( Form( "Event after calling the hadroniser (%s)", p->hadroniser->name().c_str() ) );
      if ( Logger::GetInstance()->Level>=Logger::Debug ) ev->dump();
    }
    ev->time_total = tmr.elapsed();
    p->process->addGenerationTime( ev->time_total );
    
    Debugging( Form( "Generation time:       %5.6f sec\n\t"
                     "Total time (gen+hadr): %5.6f sec",
                     ev->time_generation,
                     ev->time_total ) );

    *(p->last_event) = *( ev );
    //ev->Store(p->file);

  }

  if ( Logger::GetInstance()->Level>=Logger::DebugInsideLoop ) {
    os.str( "" ); for ( unsigned int i=0; i<ndim; i++ ) { os << Form( "%10.8f ", x[i] ); }
    Debugging( Form( "f value for dim-%d point ( %s): %4.4e", ndim, os.str().c_str(), ff ) );
  }
  
  return ff;
}

