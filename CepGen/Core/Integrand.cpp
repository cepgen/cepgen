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
                 &op2 = ev->getOneByRole( Particle::OutgoingBeam2 );
        Particles& central_system = ev->getByRole( Particle::CentralSystem );

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

        //--- prepare the function to be integrated
        p->process()->prepareKinematics();

        //--- add central system
        for ( Particles::iterator part = central_system.begin(); part != central_system.end(); ++part ) {
          part->setPdgId( p->kinematics.pair );
          part->computeMass();
        }

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
      const Particle::Momentum central_system( ev->getByRole( Particle::CentralSystem )[0].momentum() + ev->getByRole( Particle::CentralSystem )[1].momentum() );
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

