#include "CepGen/Generator.h"

namespace CepGen
{
  double
  f( double* x, size_t ndim, void* params )
  {
    std::ostringstream os;

    Parameters* p = static_cast<Parameters*>( params );
    std::shared_ptr<Event> ev = p->process()->event();

    if ( p->process()->hasEvent() ) {
      p->process()->clearEvent();

      const Particle::Momentum p1( 0., 0.,  p->kinematics.in1p ), p2( 0., 0., -p->kinematics.in2p );
      p->process()->setIncomingKinematics( p1, p2 ); // at some point introduce non head-on colliding beams?

      DebuggingInsideLoop( Form( "Function f called -- some parameters:\n\t"
                                 "  pz(p1) = %5.2f  pz(p2) = %5.2f\n\t"
                                 "  remnant mode: %d",
                                 p->kinematics.in1p, p->kinematics.in2p, p->remnant_mode ) );

      if ( p->vegas.first_run ) {

        if ( Logger::get().level >= Logger::Debug ) {
          std::ostringstream oss; oss << p->kinematics.mode;
          Debugging( Form( "Process mode considered: %s", oss.str().c_str() ) );
        }

        Particle& op1 = ev->getOneByRole( Particle::OutgoingBeam1 ),
                 &op2 = ev->getOneByRole( Particle::OutgoingBeam2 ),
                 &cs1 = ev->getOneByRole( Particle::CentralParticle1 ),
                 &cs2 = ev->getOneByRole( Particle::CentralParticle2 );

        //--- add outgoing protons or remnants
        switch ( p->kinematics.mode ) {
          case Kinematics::ElectronProton: default: { InError( Form( "Process mode %d not yet handled!", (int)p->kinematics.mode ) ); }
          case Kinematics::ElasticElastic: break; // nothing to change in the event
          case Kinematics::ElasticInelastic:
            op2.setPdgId( Particle::uQuark ); break;
          case Kinematics::InelasticElastic: // set one of the outgoing protons to be fragmented
            op1.setPdgId( Particle::uQuark ); break;
          case Kinematics::InelasticInelastic: // set both the outgoing protons to be fragmented
            op1.setPdgId( Particle::uQuark );
            op2.setPdgId( Particle::uQuark );
            break;
        }
        //PrintMessage( Form( "4 - after preparing the event kinematics: %.3e", tmr.elapsed()-now ) ); now = tmr.elapsed();

        //--- prepare the function to be integrated
        p->process()->prepareKinematics();

        //--- add central system
        cs1.setPdgId( p->kinematics.pair ); cs1.computeMass();
        cs2.setPdgId( p->kinematics.pair ); cs2.computeMass();

        p->process()->clearRun();
        p->vegas.first_run = false;
      }
    } // event is not empty

    p->process()->setPoint( ndim, x );
    if ( Logger::get().level >= Logger::DebugInsideLoop ) {
      std::ostringstream oss; for ( unsigned int i=0; i<ndim; i++ ) { oss << x[i] << " "; }
      DebuggingInsideLoop( Form( "Computing dim-%d point ( %s)", ndim, oss.str().c_str() ) );
    }

    //--- from this step on, the phase space point is supposed to be set by 'GenericProcess::setPoint()'

    p->process()->beforeComputeWeight();

    Timer tmr; // start the timer

    double integrand = p->process()->computeWeight();

    if ( integrand < 0. ) return 0.;

    //--- only fill in the process' Event object if storage is requested
    //    or if taming functions are to be applied

    if ( !p->taming_functions.empty() || p->storage() ) p->process()->fillKinematics();

    //--- once the kinematics variables have been populated,
    //    can apply the collection of taming functions

    double taming = 1.0;
    if ( p->taming_functions.has( "m_central" ) || p->taming_functions.has( "pt_central" ) ) {
      const Particle::Momentum central_system( ev->getOneByRole( Particle::CentralParticle1 ).momentum() + ev->getOneByRole( Particle::CentralParticle2 ).momentum() );
      taming *= p->taming_functions.eval( "m_central", central_system.mass() );
      taming *= p->taming_functions.eval( "pt_central", central_system.pt() );
    }
    if ( p->taming_functions.has( "q2" ) ) {
      taming *= p->taming_functions.eval( "q2", ev->getOneByRole( Particle::Parton1 ).momentum().mass() );
      taming *= p->taming_functions.eval( "q2", ev->getOneByRole( Particle::Parton2 ).momentum().mass() );
    }
    integrand *= taming;

    //--- full event content (+ hadronisation) if generating events

    if ( p->storage() ) {

      ev->time_generation = tmr.elapsed();

      //--- if an hadronisation algorithm is specified, run it

      if ( p->hadroniser() && p->kinematics.mode!=Kinematics::ElasticElastic ) {

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

      p->generation.last_event = ev;
    } // generating events

    if ( Logger::get().level>=Logger::DebugInsideLoop ) {
      os.str( "" ); for ( unsigned int i=0; i<ndim; i++ ) { os << Form( "%10.8f ", x[i] ); }
      Debugging( Form( "f value for dim-%d point ( %s): %4.4e", ndim, os.str().c_str(), integrand ) );
    }

    return integrand;
  }
}

