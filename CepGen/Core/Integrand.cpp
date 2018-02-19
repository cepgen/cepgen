#include "CepGen/Core/Timer.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Parameters.h"

#include <sstream>

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

      const Particle::Momentum p1( 0., 0.,  p->kinematics.inp.first ), p2( 0., 0., -p->kinematics.inp.second );
      p->process()->setIncomingKinematics( p1, p2 ); // at some point introduce non head-on colliding beams?

      DebuggingInsideLoop( Form( "Function f called -- some parameters:\n\t"
                                 "  pz(p1) = %5.2f  pz(p2) = %5.2f\n\t"
                                 "  remnant mode: %d",
                                 p->kinematics.inp.first, p->kinematics.inp.second, p->kinematics.structure_functions ) );

      if ( p->integrator.first_run ) {

        if ( Logger::get().level >= Logger::Debug ) {
          std::ostringstream oss; oss << p->kinematics.mode;
          Debugging( Form( "Process mode considered: %s", oss.str().c_str() ) );
        }

        Particle& op1 = ev->getOneByRole( Particle::OutgoingBeam1 ),
                 &op2 = ev->getOneByRole( Particle::OutgoingBeam2 );

        //--- add outgoing protons or remnants
        switch ( p->kinematics.mode ) {
          case Kinematics::ElectronProton: default: { InError( Form( "Process mode %d not yet handled!", (int)p->kinematics.mode ) ); }
          case Kinematics::ElasticElastic: break; // nothing to change in the event
          case Kinematics::ElasticInelastic:
            op2.setPdgId( DiffrProt, 1 ); break;
          case Kinematics::InelasticElastic: // set one of the outgoing protons to be fragmented
            op1.setPdgId( DiffrProt, 1 ); break;
          case Kinematics::InelasticInelastic: // set both the outgoing protons to be fragmented
            op1.setPdgId( DiffrProt, 1 );
            op2.setPdgId( DiffrProt, 1 );
            break;
        }

        //--- prepare the function to be integrated
        p->process()->prepareKinematics();

        //--- add central system
        Particles& central_system = ev->getByRole( Particle::CentralSystem );
        if ( central_system.size() == p->kinematics.central_system.size() ) {
          unsigned short i = 0;
          for ( Particles::iterator part = central_system.begin(); part != central_system.end(); ++part ) {
            if ( p->kinematics.central_system[i] == invalidParticle ) continue;
            part->setPdgId( p->kinematics.central_system[i] );
            part->computeMass();
            i++;
          }
        }

        p->process()->clearRun();
        p->integrator.first_run = false;
      }
    } // event is not empty

    p->process()->setPoint( ndim, x );
    if ( Logger::get().level >= Logger::DebugInsideLoop ) {
      std::ostringstream oss; for ( unsigned int i = 0; i < ndim; ++i ) { oss << x[i] << " "; }
      DebuggingInsideLoop( Form( "Computing dim-%d point ( %s)", ndim, oss.str().c_str() ) );
    }

    //--- from this step on, the phase space point is supposed to be set by 'GenericProcess::setPoint()'

    p->process()->beforeComputeWeight();

    Timer tmr; // start the timer

    double integrand = p->process()->computeWeight();

    if ( integrand <= 0. ) return 0.;

    //--- fill in the process' Event object

    p->process()->fillKinematics();

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

    if ( p->hadroniser() ) {
      double br = 0.; // branching fraction for all decays
      bool success = false;
      ev->num_hadronisation_trials = 0;
      do {
        success = p->hadroniser()->hadronise( *ev, br );
        ev->num_hadronisation_trials++;
      } while ( !success && ev->num_hadronisation_trials < p->hadroniser_max_trials );
      if ( !success ) return 0.;
      integrand *= br;
    }

    //--- apply cuts on final state system (after hadronisation!)
    //    (watch out your cuts, as this might be extremely time-consuming...)

    if ( p->kinematics.cuts.central_particles.size() > 0 ) {
      std::map<ParticleCode,std::map<Cuts::Central,Kinematics::Limits> >::const_iterator it_c;
      const Particles cs = ev->getByRole( Particle::CentralSystem );
      for ( Particles::const_iterator it_p = cs.begin(); it_p != cs.end(); ++it_p ) {
        it_c = p->kinematics.cuts.central_particles.find( it_p->pdgId() );
        if ( it_c == p->kinematics.cuts.central_particles.end() ) continue;
        // retrieve all cuts associated to this final state particle
        const std::map<Cuts::Central,Kinematics::Limits>& cm = it_c->second;
        // apply these cuts on the given particle
        if ( cm.count( Cuts::pt_single ) > 0 && !cm.at( Cuts::pt_single ).passes( it_p->momentum().pt() ) ) return 0.;
        //std::cout << it_c->first << "\t" << it_p->momentum().pt() << "\t" << cm.at( Cuts::pt_single ).passes( it_p->momentum().pt() ) << std::endl;
        if ( cm.count( Cuts::energy_single ) > 0 && !cm.at( Cuts::energy_single ).passes( it_p->momentum().energy() ) ) return 0.;
        if ( cm.count( Cuts::eta_single ) > 0 && !cm.at( Cuts::eta_single ).passes( it_p->momentum().eta() ) ) return 0.;
        if ( cm.count( Cuts::rapidity_single ) > 0 && !cm.at( Cuts::rapidity_single ).passes( it_p->momentum().rapidity() ) ) return 0.;
      }
    }

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

    if ( Logger::get().level >= Logger::DebugInsideLoop ) {
      os.str( "" ); for ( unsigned int i = 0; i < ndim; ++i ) { os << Form( "%10.8f ", x[i] ); }
      Debugging( Form( "f value for dim-%d point ( %s): %4.4e", ndim, os.str().c_str(), integrand ) );
    }

    return integrand;
  }
}

