#include "CepGenAddOns/EventInterfaces/HepMCEventInterface.h"

#include "CepGen/Parameters.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"

#ifdef HEPMC3
# include "HepMC3/Version.h"
# include "HepMC3/FourVector.h"
# include "HepMC3/GenEvent.h"
# include "HepMC3/GenRunInfo.h"
# include "HepMC3/GenVertex.h"
# include "HepMC3/GenParticle.h"
#else
# include "HepMC/Version.h"
# if !defined( HEPMC_VERSION_CODE ) // HepMC v2
#  include "HepMC/SimpleVector.h"
#  include "HepMC/GenEvent.h"
#  include "HepMC/GenVertex.h"
#  include "HepMC/GenParticle.h"
#  define BUILD( type ) new type
# else
#  include "HepMC/FourVector.h"
#  include "HepMC/GenEvent.h"
#  include "HepMC/GenRunInfo.h"
#  include "HepMC/GenVertex.h"
#  include "HepMC/GenParticle.h"
#  define HEPMC3
# endif
#endif
#ifndef BUILD
# define BUILD( type ) make_shared<type>
#endif

#include <list>
#include <numeric>

namespace HepMC
{
  CepGenEvent::CepGenEvent( const cepgen::Event& evt ) :
    GenEvent( Units::GEV, Units::MM )
  {
#ifdef HEPMC3
    add_attribute( "AlphaQCD", make_shared<DoubleAttribute>( cepgen::constants::ALPHA_QCD ) );
    add_attribute( "AlphaEM", make_shared<DoubleAttribute>( cepgen::constants::ALPHA_EM ) );
#else
    set_alphaQCD( cepgen::constants::ALPHA_QCD );
    set_alphaQED( cepgen::constants::ALPHA_EM );
#endif

    weights().push_back( 1. ); // unweighted events

    // filling the particles content
    const FourVector origin( 0., 0., 0., 0. );
    int cm_id = 0;

    auto v1 = BUILD(GenVertex)( origin ), v2 = BUILD(GenVertex)( origin ), vcm = BUILD(GenVertex)( origin );
    unsigned short idx = 0;
    for ( const auto& part_orig : evt.particles() ) {
      const auto& mom_orig = part_orig.momentum();
      FourVector pmom( mom_orig.px(), mom_orig.py(), mom_orig.pz(), part_orig.energy() );
      auto part = BUILD(GenParticle)( pmom, part_orig.integerPdgId(), (int)part_orig.status() );
#ifndef HEPMC3
      part->suggest_barcode( idx );
#endif
      part->set_generated_mass( cepgen::PDG::get().mass( part_orig.pdgId() ) );
      assoc_map_[idx] = part;

      switch ( part_orig.role() ) {
        case cepgen::Particle::IncomingBeam1:
          v1->add_particle_in( part );
          break;
        case cepgen::Particle::IncomingBeam2:
          v2->add_particle_in( part );
          break;
        case cepgen::Particle::OutgoingBeam1:
          v1->add_particle_out( part );
          break;
        case cepgen::Particle::OutgoingBeam2:
          v2->add_particle_out( part );
          break;
        case cepgen::Particle::Parton1:
          v1->add_particle_out( part );
          vcm->add_particle_in( part );
          break;
        case cepgen::Particle::Parton2:
          v2->add_particle_out( part );
          vcm->add_particle_in( part );
          break;
        case cepgen::Particle::Intermediate:
          // skip the two-parton system and propagate the parentage
          cm_id = idx;
          continue;
        case cepgen::Particle::CentralSystem: default: {
          const auto& moth = part_orig.mothers();
          if ( moth.empty() )
            // skip disconnected lines
            continue;
          // get mother(s) id(s)
          const short m1 = *moth.begin();
          const short m2 = moth.size() > 1 ? *moth.rbegin() : -1;
          // check if particle is connected to the two-parton system
          if ( m1 == cm_id
           || ( m2 >= 0 && ( m1 < cm_id && cm_id <= m2 ) ) ) // also supports range
            vcm->add_particle_out( part );
          // if part of the decay chain of central system, find parents
          else if ( assoc_map_.count( m1 ) != 0 ) {
            auto vprod = assoc_map_.at( m1 )->end_vertex();
            std::list<short> ids{ m1 }; // list of mother particles
            if ( assoc_map_.count( m2 ) != 0 && m2 > m1 ) {
              ids.resize( m2-m1+1 );
              std::iota( ids.begin(), ids.end(), m1 );
            }
            if ( !vprod ) {
              vprod = BUILD(GenVertex)();
              for ( const auto& id : ids )
                vprod->add_particle_in( assoc_map_.at( id ) );
              add_vertex( vprod );
            }
            vprod->add_particle_out( part );
          }
          else
            throw CG_FATAL( "HepMCHandler:fillEvent" )
              << "Other particle requested! Not yet implemented!";
        } break;
      }
      idx++;
    }
    add_vertex( v1 );
    add_vertex( v2 );
    add_vertex( vcm );

#ifndef HEPMC3
    set_beam_particles( *v1->particles_in_const_begin(), *v2->particles_in_const_begin() );
    set_signal_process_vertex( vcm );
#endif
  }
}

#undef BUILD
