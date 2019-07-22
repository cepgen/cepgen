#include "CepGen/IO/HepMCEventInterface.h"

#include "CepGen/Parameters.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"

#ifdef HEPMC3
#  include "HepMC3/Version.h"
#  include "HepMC3/FourVector.h"
#  include "HepMC3/GenEvent.h"
#  include "HepMC3/GenRunInfo.h"
#  include "HepMC3/GenVertex.h"
#  include "HepMC3/GenParticle.h"
namespace HepMC3
{
#else
#  include "HepMC/Version.h"
#  if !defined( HEPMC_VERSION_CODE ) // HepMC v2
#    include "HepMC/SimpleVector.h"
#    include "HepMC/GenEvent.h"
#    include "HepMC/GenVertex.h"
#    include "HepMC/GenParticle.h"
#  else
#    include "HepMC/FourVector.h"
#    include "HepMC/GenEvent.h"
#    include "HepMC/GenRunInfo.h"
#    include "HepMC/GenVertex.h"
#    include "HepMC/GenParticle.h"
#    define HEPMC3
#  endif
namespace HepMC
{
#endif

  CepGenEvent::CepGenEvent() :
    GenEvent( Units::GEV, Units::MM )
  {
#ifdef HEPMC3
    GenEvent::add_attribute( "AlphaQCD", make_shared<DoubleAttribute>( cepgen::constants::ALPHA_QCD ) );
    GenEvent::add_attribute( "AlphaEM", make_shared<DoubleAttribute>( cepgen::constants::ALPHA_EM ) );
#else
    GenEvent::set_alphaQCD( cepgen::constants::ALPHA_QCD );
    GenEvent::set_alphaQED( cepgen::constants::ALPHA_EM );
#endif
  }

  void
  CepGenEvent::feedEvent( const cepgen::Event& evt )
  {
    GenEvent::clear();
    GenEvent::weights().push_back( 1. ); // unweighted events

    // filling the particles content
    const FourVector origin( 0., 0., 0., 0. );
    auto& part_vec = evt.particles();

    int cm_id = 0, idx = 1;

#ifdef HEPMC3
    GenVertexPtr v1 = make_shared<GenVertex>( origin ),
                 v2 = make_shared<GenVertex>( origin ),
                 vcm = make_shared<GenVertex>( origin );
#else
    GenVertex* v1 = new GenVertex( origin ),
              *v2 = new GenVertex( origin ),
              *vcm = new GenVertex( origin );
#endif
    for ( unsigned int i = 0; i < part_vec.size(); ++i ) {
      const auto& part_orig = part_vec.at( i );
      FourVector pmom( part_orig.momentum().px(),
                       part_orig.momentum().py(),
                       part_orig.momentum().pz(),
                       part_orig.energy() );
#ifdef HEPMC3
      auto part = make_shared<GenParticle>( pmom, part_orig.integerPdgId(), (int)part_orig.status() );
#else
      auto part = new GenParticle( pmom, part_orig.integerPdgId(), (int)part_orig.status() );
      part->suggest_barcode( idx++ );
#endif
      const auto& moth = part_orig.mothers();

      switch ( part_orig.role() ) {
        case cepgen::Particle::IncomingBeam1: { v1->add_particle_in( part ); } break;
        case cepgen::Particle::IncomingBeam2: { v2->add_particle_in( part ); } break;
        case cepgen::Particle::OutgoingBeam1: { v1->add_particle_out( part ); } break;
        case cepgen::Particle::OutgoingBeam2: { v2->add_particle_out( part ); } break;
        case cepgen::Particle::Parton1:       { v1->add_particle_out( part ); vcm->add_particle_in( part ); } break;
        case cepgen::Particle::Parton2:       { v2->add_particle_out( part ); vcm->add_particle_in( part ); } break;
        case cepgen::Particle::Intermediate:  { cm_id = i; continue; } break;
        case cepgen::Particle::CentralSystem:
        default: {
          if ( moth.empty() ) continue;
          if ( *moth.begin() == cm_id ) vcm->add_particle_out( part );
          else
            throw CG_FATAL( "HepMCHandler:fillEvent" )
              << "Other particle requested! Not yet implemented!";
        } break;
      }
      idx++;
    }
    GenEvent::add_vertex( v1 );
    GenEvent::add_vertex( v2 );
    GenEvent::add_vertex( vcm );

#ifndef HEPMC3
    GenEvent::set_beam_particles( *v1->particles_in_const_begin(), *v2->particles_in_const_begin() );
    GenEvent::set_signal_process_vertex( *v1->vertices_begin() );
    GenEvent::set_beam_particles( *v1->particles_in_const_begin(), *v2->particles_in_const_end() );
#endif
  }
}

