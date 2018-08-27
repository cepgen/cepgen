#include "CepGen/Core/Timer.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/TamingFunction.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#include "CepGen/Parameters.h"

#include <sstream>
#include <fstream>

namespace CepGen
{
  namespace Integrand
  {
    Logger::Level log_level;
    Timer tmr;

    double
    eval( double* x, size_t ndim, void* params )
    {
      log_level = Logger::get().level;
      std::shared_ptr<Event> ev;

      Parameters* p = static_cast<Parameters*>( params );
      if ( !p )
        throw CG_FATAL( "Integrand" ) << "Failed to retrieve the run parameters!";

      Process::GenericProcess* proc = p->process();
      if ( !proc )
        throw CG_FATAL( "Integrand" ) << "Failed to retrieve the process!";

      //=============================================================================================
      // start the timer
      //=============================================================================================

      tmr.reset();

      //=============================================================================================
      // prepare the event content prior to the process generation
      //=============================================================================================

      if ( proc->hasEvent() ) {
        ev = proc->event();

        if ( proc->first_run ) {
          CG_DEBUG( "Integrand" )
            << "Computation launched for " << p->processName() << " process "
            << "0x" << std::hex << p->process() << std::dec << ".\n\t"
            << "Process mode considered: " << p->kinematics.mode << "\n\t"
            << "  pz(p1) = " << p->kinematics.incoming_beams.first.pz << "\n\t"
            << "  pz(p2) = " << p->kinematics.incoming_beams.second.pz << "\n\t"
            << "  structure functions: " << p->kinematics.structure_functions;

          p->clearRunStatistics();
          proc->first_run = false;
        } // passed the first-run preparation

        proc->clearEvent();
      } // event is not empty

      //=============================================================================================
      // specify the phase space point to probe
      //=============================================================================================

      proc->setPoint( ndim, x );

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
          if ( !cuts_pdgid.pt_single.passes( part.momentum().pt() ) )
            return 0.;
          if ( !cuts_pdgid.energy_single.passes( part.momentum().energy() ) )
            return 0.;
          if ( !cuts_pdgid.eta_single.passes( part.momentum().eta() ) )
            return 0.;
          if ( !cuts_pdgid.rapidity_single.passes( part.momentum().rapidity() ) )
            return 0.;
        }
      }

      //=============================================================================================
      // store the last event in the parameters for its usage by the end user
      //=============================================================================================

      if ( p->storage() ) {
        p->process()->last_event = ev;
        p->process()->last_event->time_total = tmr.elapsed();

        CG_DEBUG( "Integrand" )
          << "[process 0x" << std::hex << p->process() << std::dec << "] "
          << "Individual time (gen+hadr+cuts): " << p->process()->last_event->time_total*1.e3 << " ms";
      }

      //=============================================================================================
      // a bit of useful debugging
      //=============================================================================================

      if ( CG_EXCEPT_MATCH( "Integrand", debugInsideLoop ) ) {
        std::ostringstream oss;
        for ( unsigned short i = 0; i < ndim; ++i )
          oss << Form( "%10.8f ", x[i] );
        CG_DEBUG( "Integrand" )
          << "f value for dim-" << ndim << " point ( " << oss.str() << "): "
          << integrand;
      }

      return integrand;
    }
  }
}

