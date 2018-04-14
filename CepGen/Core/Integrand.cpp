#include "CepGen/Core/Timer.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/TamingFunction.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/Kinematics.h"

#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#include "CepGen/Parameters.h"

#include <sstream>
#include <fstream>

namespace CepGen
{
  namespace Integrand
  {
    Logger::LoggingLevel log_level;
    Timer tmr;

    double
    eval( double* x, size_t ndim, void* params )
    {
      log_level = Logger::get().level;
      std::shared_ptr<Event> ev;

      Parameters* p = static_cast<Parameters*>( params );
      if ( !p )
        throw Exception( __PRETTY_FUNCTION__, "Failed to retrieve the run parameters!", FatalError );

      Process::GenericProcess* proc = p->process();
      if ( !proc )
        throw Exception( __PRETTY_FUNCTION__, "Failed to retrieve the process!", FatalError );

      //=============================================================================================
      // start the timer
      //=============================================================================================

      tmr.reset();

      //=============================================================================================
      // prepare the event content prior to the process generation
      //=============================================================================================

      if ( proc->hasEvent() ) {
        ev = proc->event();
        proc->clearEvent();

        if ( proc->first_run ) {
          if ( log_level >= Logger::Debug )
            Debugging( Form( "Computation launched for %s process 0x%zx", p->processName().c_str(), p->process() ) );

          const Particle::Momentum p1( 0., 0.,  p->kinematics.inp.first ), p2( 0., 0., -p->kinematics.inp.second );
          proc->setIncomingKinematics( p1, p2 ); // at some point introduce non head-on colliding beams?

          if ( log_level >= Logger::Debug ) {
            std::ostringstream oss; oss << p->kinematics.mode;
            Debugging( Form( "Process mode considered: %s\n\t"
                             "  pz(p1) = %5.2f\n\t"
                             "  pz(p2) = %5.2f\n\t"
                             "  remnant mode: %d",
                              oss.str().c_str(),
                              p->kinematics.inp.first, p->kinematics.inp.second,
                              p->kinematics.structure_functions ) );
          }

          //=========================================================================================
          // prepare the function to be integrated
          //=========================================================================================

          proc->prepareKinematics();

          //=========================================================================================
          // add central system
          //=========================================================================================

          Particles& central_system = ev->getByRole( Particle::CentralSystem );
          if ( central_system.size() == p->kinematics.central_system.size() ) {
            unsigned short i = 0;
            for ( auto& part : central_system ) {
              if ( p->kinematics.central_system[i] == invalidParticle )
                continue;
              part.setPdgId( p->kinematics.central_system[i] );
              part.computeMass();
              i++;
            }
          }

          p->clearRunStatistics();
          proc->first_run = false;
        } // passed the first-run preparation
      } // event is not empty

      //=============================================================================================
      // specify the phase space point to probe
      //=============================================================================================

      proc->setPoint( ndim, x );
      if ( log_level >= Logger::DebugInsideLoop ) {
        std::ostringstream oss;
        for ( unsigned int i = 0; i < ndim; ++i )
          oss << ( i == 0 ? "" : " " ) << x[i];
        DebuggingInsideLoop( Form( "Computing dim-%d point (%s)", ndim, oss.str().c_str() ) );
      }

      //=============================================================================================
      // from this step on, the phase space point is supposed to be set
      //=============================================================================================

      p->process()->beforeComputeWeight();

      double integrand = p->process()->computeWeight();

      //=============================================================================================
      // invalidate any unphysical behaviour
      //=============================================================================================

      if ( integrand <= 0. )
        return 0.;

      //=============================================================================================
      // speed up the integration process if no event needs to be generated
      //=============================================================================================

      if ( !p->storage()
        && !p->taming_functions
        && !p->hadroniser()
        &&  p->kinematics.cuts.central_particles.size() == 0 )
        return integrand;
//std::cout << integrand << std::endl;

      //=============================================================================================
      // fill in the process' Event object
      //=============================================================================================

      p->process()->fillKinematics();

      //=============================================================================================
      // once the kinematics variables have been populated, can apply the collection of taming functions
      //=============================================================================================

      if ( p->taming_functions ) {
        if ( p->taming_functions->has( "m_central" )
          || p->taming_functions->has( "pt_central" ) ) {

          // build the kinematics of the central system
          Particle::Momentum central_system;
          for ( const auto& part : ev->getByRole( Particle::CentralSystem ) )
            central_system += part.momentum();

          // tame the cross-section by the reweighting function
          if ( p->taming_functions->has( "m_central" ) )
            integrand *= p->taming_functions->eval( "m_central", central_system.mass() );
          if ( p->taming_functions->has( "pt_central" ) )
            integrand *= p->taming_functions->eval( "pt_central", central_system.pt() );
        }
        if ( p->taming_functions->has( "q2" ) ) {
          integrand *= p->taming_functions->eval( "q2", -ev->getOneByRole( Particle::Parton1 ).momentum().mass() );
          integrand *= p->taming_functions->eval( "q2", -ev->getOneByRole( Particle::Parton2 ).momentum().mass() );
        }
      }

      if ( integrand <= 0. )
        return 0.;

      //=============================================================================================
      // set the CepGen part of the event generation
      //=============================================================================================

      if ( p->storage() )
        ev->time_generation = tmr.elapsed();

      //=============================================================================================
      // event hadronisation and resonances decay
      //=============================================================================================

      if ( p->hadroniser() ) {
        double br = -1.;
        if ( !p->hadroniser()->run( *ev, br, p->storage() ) || br == 0. )
//        if ( !p->hadroniser()->run( *ev, br, true ) || br == 0. )
          return 0.;
        integrand *= br; // branching fraction for all decays
      }

      //=============================================================================================
      // apply cuts on final state system (after hadronisation!)
      // (watch out your cuts, as this might be extremely time-consuming...)
      //=============================================================================================

      if ( p->kinematics.cuts.central_particles.size() > 0 ) {
        for ( const auto& part : ev->getByRole( Particle::CentralSystem ) ) {
          // retrieve all cuts associated to this final state particle
          if ( p->kinematics.cuts.central_particles.count( part.pdgId() ) == 0 )
            continue;
          const auto& cuts_pdgid = p->kinematics.cuts.central_particles.at( part.pdgId() );
          // apply these cuts on the given particle
          if ( cuts_pdgid.count( Cuts::pt_single ) > 0
           && !cuts_pdgid.at( Cuts::pt_single ).passes( part.momentum().pt() ) )
            return 0.;
          if ( cuts_pdgid.count( Cuts::energy_single ) > 0
           && !cuts_pdgid.at( Cuts::energy_single ).passes( part.momentum().energy() ) )
            return 0.;
          if ( cuts_pdgid.count( Cuts::eta_single ) > 0
           && !cuts_pdgid.at( Cuts::eta_single ).passes( part.momentum().eta() ) )
            return 0.;
          if ( cuts_pdgid.count( Cuts::rapidity_single ) > 0
           && !cuts_pdgid.at( Cuts::rapidity_single ).passes( part.momentum().rapidity() ) )
            return 0.;
        }
      }

      //=============================================================================================
      // store the last event in the parameters for its usage by the end user
      //=============================================================================================

      if ( p->storage() ) {
        p->process()->last_event = ev;
        p->process()->last_event->time_total = tmr.elapsed();

        if ( log_level >= Logger::Debug )
          Debugging( Form( "[process 0x%zx] Individual time (gen+hadr+cuts): %5.6f ms",
                           p->process(), p->process()->last_event->time_total*1.e3 ) );
      }

      //=============================================================================================
      // a bit of useful debugging
      //=============================================================================================

      if ( log_level >= Logger::DebugInsideLoop ) {
        std::ostringstream oss;
        for ( unsigned short i = 0; i < ndim; ++i )
          oss << Form( "%10.8f ", x[i] );
        Debugging( Form( "f value for dim-%d point ( %s): %4.4e",
                         ndim, oss.str().c_str(), integrand ) );
      }

      return integrand;
    }
  }
}

