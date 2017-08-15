#include "CepGen/Generator.h"

namespace CepGen
{
  double
  f( double* x, size_t ndim, void* params )
  {
    Timer tmr;
    std::ostringstream os;

    Parameters* p = static_cast<Parameters*>( params );

    //float now = tmr.elapsed();
    const Particle::Momentum p1( 0., 0.,  p->kinematics.in1p ), p2( 0., 0., -p->kinematics.in2p );
std::cout << p1 << "\t" << p2 << std::endl;
    p->process()->setIncomingKinematics( p1, p2 ); // at some point introduce non head-on colliding beams?
    //PrintMessage( Form( "0 - after setting the kinematics: %.3e", tmr.elapsed()-now ) ); now = tmr.elapsed();

    DebuggingInsideLoop( Form( "Function f called -- some parameters:\n\t"
                               "  pz(p1) = %5.2f  pz(p2) = %5.2f\n\t"
                               "  remnant mode: %d",
                               p->kinematics.in1p, p->kinematics.in2p, p->remnant_mode ) );

    p->process()->clearEvent();
    //PrintMessage( Form( "1 - after clearing the event: %.3e", tmr.elapsed()-now ) ); now = tmr.elapsed();
  
    std::shared_ptr<Event> ev = p->process()->event();

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

    if ( p->generation ) { // MC events generation
      p->process()->fillKinematics( false );

      ev->time_generation = tmr.elapsed();

      if ( p->hadroniser() && p->process_mode!=Kinematics::ElasticElastic ) {

        Debugging( Form( "Event before calling the hadroniser (%s)", p->hadroniser()->name().c_str() ) );
        if ( Logger::get().level>=Logger::Debug ) ev->dump();

        unsigned int num_hadr_trials = 0;
        bool hadronised = false;
        while ( !hadronised && num_hadr_trials <= p->hadroniser_max_trials ) {
          try {
            hadronised = p->hadroniser()->hadronise( ev.get() );
          } catch ( Exception& e ) { e.dump(); }

          if ( num_hadr_trials>0 ) { Debugging( Form( "Hadronisation failed. Trying for the %dth time", num_hadr_trials+1 ) ); }
          num_hadr_trials++;
        }
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

      *( p->last_event ) = *( ev );
    }

    if ( Logger::get().level>=Logger::DebugInsideLoop ) {
      os.str( "" ); for ( unsigned int i=0; i<ndim; i++ ) { os << Form( "%10.8f ", x[i] ); }
      Debugging( Form( "f value for dim-%d point ( %s): %4.4e", ndim, os.str().c_str(), ff ) );
    }

    return ff;
  }
}

