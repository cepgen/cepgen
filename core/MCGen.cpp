#include "MCGen.h"

MCGen::MCGen() :
  fVegas( 0 ), fCrossSection( -1. ), fCrossSectionError( -1. )
{
  Debugging( "Generator initialized" );
  
  try { PrintHeader(); } catch ( Exception& e ) { e.Dump(); }
  
  srand( time( 0 ) ); // Random number initialization
  
  this->parameters = new Parameters;
}

MCGen::MCGen( Parameters *ip_ ) :
  parameters( ip_ ), fVegas( 0 )
{}

MCGen::~MCGen()
{
  if ( parameters->generation and parameters->process and parameters->process->NumGeneratedEvents()>0 ) {
    Information( Form( "Mean generation time / event: %.3f ms", parameters->process->TotalGenerationTime()*1.e3/parameters->process->NumGeneratedEvents() ) );
  }

  if ( fVegas ) delete fVegas;
  if ( parameters ) delete parameters;
}

void
MCGen::PrintHeader()
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
MCGen::BuildVegas()
{
  if ( Logger::GetInstance()->Level>=Logger::Debug ) {
    std::ostringstream topo; topo << parameters->process_mode;
    Debugging( Form( "Considered topology: %s case", topo.str().c_str() ) );
  }
  
  if ( fVegas ) delete fVegas;
  fVegas = new Vegas( GetNdim(), f, parameters );
}

void
MCGen::ComputeXsection( double* xsec_, double *err_ )
{
  if ( !fVegas ) BuildVegas();

  Information( "Starting the computation of the process cross-section" );

  try { PrepareFunction(); } catch ( Exception& e ) { e.Dump(); }
  fVegas->Integrate( xsec_, err_ );
  
  fCrossSection = *xsec_;
  fCrossSectionError = *err_;
  fHasCrossSection = true;
  
  Information( Form( "Total cross section: %f +/- %f pb", *xsec_, *err_ ) );
}

Event*
MCGen::GenerateOneEvent()
{
  bool good = false;
  if ( !fHasCrossSection ) {
    double xsec, err;
    ComputeXsection( &xsec, &err );
  }
  while ( !good ) { good = fVegas->GenerateOneEvent(); }

  last_event = this->parameters->last_event;
  return static_cast<Event*>( last_event );
}

void
MCGen::PrepareFunction()
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
  kin.emin = parameters->minenergy;
  kin.emax = parameters->maxenergy;
  kin.mxmin = parameters->minmx;
  kin.mxmax = parameters->maxmx;
  parameters->process->AddEventContent();
  parameters->process->SetKinematics( kin );
  Debugging( "Function prepared to be integrated!" );
}

double f( double* x_, size_t ndim_, void* params_ )
{
  double ff;
  Parameters *p;
  Timer tmr;
  bool hadronised;
  double num_hadr_trials;
  std::ostringstream os;

  p = static_cast<Parameters*>( params_ );

  //FIXME at some point introduce non head-on colliding beams ?

  const Particle::Momentum p1( 0., 0.,  p->in1p ),
                           p2( 0., 0., -p->in2p );
  p->process->SetIncomingKinematics( p1, p2 );
  p->process->SetPoint( ndim_, x_ );

  if ( Logger::GetInstance()->Level>=Logger::DebugInsideLoop ) {
    os.str(""); for ( unsigned int i=0; i<ndim_; i++ ) { os << x_[i] << " "; }
    DebuggingInsideLoop( Form( "Computing dim-%d point ( %s)", ndim_, os.str().c_str() ) );
  }

  tmr.reset();

  ff = 0.;

  DebuggingInsideLoop( Form( "Function f called -- some parameters:\n\t"
                             "  pz(p1) = %5.2f  pz(p2) = %5.2f\n\t"
                             "  remnant mode: %d",
                             p->in1p, p->in2p, p->remnant_mode ) );
    
  //float now = tmr.elapsed();
  //std::cout << "0: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
  p->process->ClearEvent();
  //std::cout << "1: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
  
  Event* ev = p->process->GetEvent();
  
  //std::cout << "2: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();*/

  if ( p->first_run ) {
    // Then add outgoing leptons
    ev->GetOneByRole( Particle::CentralParticle1 )->SetPDGId( p->pair );
    ev->GetOneByRole( Particle::CentralParticle2 )->SetPDGId( p->pair );

    //std::cout << "3: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
    // Then add outgoing protons or remnants
    switch ( p->process_mode ) {
      case Kinematics::ElectronProton: { InError( "Not handled yet!" ); }
      case Kinematics::ElasticElastic: break; // nothing to change in the event
      case Kinematics::ElasticInelastic:
      case Kinematics::InelasticElastic: // set one of the outgoing protons to be fragmented
        ev->GetOneByRole( Particle::OutgoingBeam1 )->SetPDGId( Particle::uQuark ); break;
      case Kinematics::InelasticInelastic: // set both the outgoing protons to be fragmented
        ev->GetOneByRole( Particle::OutgoingBeam1 )->SetPDGId( Particle::uQuark );
        ev->GetOneByRole( Particle::OutgoingBeam2 )->SetPDGId( Particle::uQuark );
        break;
    }
    
    //std::cout << "4: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
    // Prepare the function to be integrated
    p->process->PrepareKinematics();
    p->process->ClearRun();
    p->first_run = false;
  }
   
  p->process->BeforeComputeWeight();

  //std::cout << "5: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
  ff = p->process->ComputeWeight();
  //std::cout << "6: " << (tmr.elapsed()-now) << std::endl; now = tmr.elapsed();
  if (ff<0.) return 0.;
  
  if (p->store) { // MC events generation
    p->process->FillKinematics( false );
    
    ev->time_generation = tmr.elapsed();

    if ( p->hadroniser and p->process_mode!=Kinematics::ElasticElastic ) {
      
      Debugging( Form( "Event before calling the hadroniser (%s)", p->hadroniser->GetName().c_str() ) );
      if ( Logger::GetInstance()->Level>=Logger::Debug ) ev->Dump();
      
      num_hadr_trials = 0;
      do {
        try { hadronised = p->hadroniser->Hadronise( ev ); } catch ( Exception& e ) { e.Dump(); }

        if ( num_hadr_trials>0 ) { Debugging( Form( "Hadronisation failed. Trying for the %dth time", num_hadr_trials+1 ) ); }
        
        num_hadr_trials++;
      } while ( !hadronised and num_hadr_trials<=p->hadroniser_max_trials );
      if ( !hadronised ) return 0.; //FIXME
      
      ev->num_hadronisation_trials = num_hadr_trials;

      Debugging( Form( "Event hadronisation succeeded after %d trial(s)", ev->num_hadronisation_trials ) );

      if ( num_hadr_trials>p->hadroniser_max_trials ) return 0.; //FIXME
      
      Debugging( Form( "Event after calling the hadroniser (%s)", p->hadroniser->GetName().c_str() ) );
      if ( Logger::GetInstance()->Level>=Logger::Debug ) ev->Dump();
    }
    ev->time_total = tmr.elapsed();
    p->process->AddGenerationTime( ev->time_total );
    
    Debugging( Form( "Generation time:       %5.6f sec\n\t"
                     "Total time (gen+hadr): %5.6f sec",
                     ev->time_generation,
                     ev->time_total ) );

    *(p->last_event) = *( ev );
    //ev->Store(p->file);

  }

  if ( Logger::GetInstance()->Level>=Logger::DebugInsideLoop ) {
    os.str( "" ); for ( unsigned int i=0; i<ndim_; i++ ) { os << Form( "%10.8f ", x_[i] ); }
    Debugging( Form( "f value for dim-%d point ( %s): %4.4e", ndim_, os.str().c_str(), ff ) );
  }
  
  return ff;
}

