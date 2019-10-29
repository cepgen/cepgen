#include "CepGen/Core/Timer.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/EventBrowser.h"

#include "CepGen/Physics/TamingFunction.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/GenericProcess.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/GenericExportHandler.h"

#include "CepGen/Parameters.h"

#include <cmath>
#include <sstream>
#include <fstream>

namespace cepgen
{
  namespace integrand
  {
    utils::Logger::Level log_level;
    utils::Timer tmr;

    double
    eval( double* x, size_t ndim, void* func_params )
    {
      log_level = utils::Logger::get().level;
      std::shared_ptr<Event> ev;

      Parameters* params = nullptr;
      if ( !func_params || !( params = static_cast<Parameters*>( func_params ) ) )
        throw CG_FATAL( "Integrand" ) << "Failed to retrieve the run parameters!";

      proc::GenericProcess* proc = params->process();
      if ( !proc )
        throw CG_FATAL( "Integrand" ) << "Failed to retrieve the process!";

      //================================================================
      // start the timer
      //================================================================

      tmr.reset();

      //================================================================
      // prepare the event content prior to the process generation
      //================================================================

      if ( proc->hasEvent() ) // event is not empty
        ev = proc->event();

      params->prepareRun();

      //================================================================
      // specify the phase space point to probe
      //================================================================

      proc->setPoint( ndim, x );

      //================================================================
      // from this step on, the phase space point is supposed to be set
      //================================================================

      proc->beforeComputeWeight();
      double weight = proc->computeWeight();

      //================================================================
      // invalidate any unphysical behaviour
      //================================================================

      if ( weight <= 0. )
        return 0.;

      //================================================================
      // speed up the integration process if no event is to be generated
      //================================================================

      if ( !ev )
        return weight;

      if ( !params->storage()
        && !params->taming_functions.empty()
        && !params->eventModifiersSequence().empty()
        &&  params->kinematics.cuts.central_particles.empty() )
        return weight;

      //================================================================
      // fill in the process' Event object
      //================================================================

      proc->fillKinematics();

      //================================================================
      // once the kinematics variables have been populated, can apply
      // the collection of taming functions
      //================================================================

      try {
        utils::EventBrowser bws;
        for ( const auto& tam : params->taming_functions )
          weight *= tam.function.eval( bws.get( *ev, tam.var_orig ) );
      } catch ( const Exception& ) {
        throw CG_FATAL( "Integrand" )
          << "Failed to apply taming function(s) taming!";
      }

      if ( weight <= 0. )
        return 0.;

      //================================================================
      // set the CepGen part of the event generation
      //================================================================

      if ( params->storage() )
        ev->time_generation = tmr.elapsed();

      //================================================================
      // event hadronisation and resonances decay
      //================================================================

      if ( !params->eventModifiersSequence().empty() ) {
        double br = -1.;
        for ( auto& mod : params->eventModifiersSequence() ) {
          if ( !mod->run( *ev, br, params->storage() ) || br == 0. )
            return 0.;
          weight *= br; // branching fraction for all decays
        }
      }

      //================================================================
      // apply cuts on final state system (after hadronisation!)
      // (polish your cuts, as this might be very time-consuming...)
      //================================================================

      if ( !params->kinematics.cuts.central_particles.empty() ) {
        for ( const auto& part : (*ev)[Particle::CentralSystem] ) {
          // retrieve all cuts associated to this final state particle
          if ( params->kinematics.cuts.central_particles.count( part.pdgId() ) == 0 )
            continue;
          const auto& cuts_pdgid = params->kinematics.cuts.central_particles.at( part.pdgId() );
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
      for ( const auto& system : { Particle::OutgoingBeam1, Particle::OutgoingBeam2 } )
        for ( const auto& part : (*ev)[system] ) {
          if ( part.status() != Particle::Status::FinalState )
            continue;
          if ( !params->kinematics.cuts.remnants.energy_single.passes( part.momentum().energy() ) )
            return 0.;
          if ( !params->kinematics.cuts.remnants.rapidity_single.passes( fabs( part.momentum().rapidity() ) ) )
            return 0.;
        }

      //================================================================
      // store the last event in parameters block for a later usage
      //================================================================

      if ( params->storage() ) {
        ev->weight = weight;
        proc->last_event = ev;
        proc->last_event->time_total = tmr.elapsed();

        CG_DEBUG( "Integrand" )
          << "[process 0x" << std::hex << proc << std::dec << "] "
          << "Individual time (gen+hadr+cuts): " << proc->last_event->time_total*1.e3 << " ms";
      }

      //================================================================
      // a bit of useful debugging
      //================================================================

      if ( CG_LOG_MATCH( "Integrand", debugInsideLoop ) ) {
        std::ostringstream oss;
        for ( unsigned short i = 0; i < ndim; ++i )
          oss << Form( "%10.8f ", x[i] );
        CG_DEBUG( "Integrand" )
          << "f value for dim-" << ndim << " point ( " << oss.str() << "): "
          << weight;
      }

      return weight;
    }
  }
}
