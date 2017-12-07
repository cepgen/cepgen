#include "HepMCHandler.h"
#include "CepGen/Parameters.h"
#include "CepGen/Core/Exception.h"

#ifdef LIBHEPMC

#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"

using namespace CepGen::OutputHandler;

HepMCHandler::HepMCHandler( const char* filename, const ExportHandler::OutputType& type ) :
  ExportHandler( type ),
  event( new HepMC::GenEvent() )
{
#ifdef HEPMC_VERSION3
  output = std::unique_ptr<HepMC::WriterAscii>( new HepMC::WriterAscii( filename ) );
#else
  output = std::unique_ptr<HepMC::IO_GenEvent>( new HepMC::IO_GenEvent( filename ) );
#endif
}

HepMCHandler::~HepMCHandler()
{}

void
HepMCHandler::operator<<( const Event* evt )
{
  fillEvent( evt );
  if ( !event.get() ) {
    throw Exception( __PRETTY_FUNCTION__, "Failed to retrieve the HepMC event to be stored!", FatalError );
  }
#ifdef HEPMC_VERSION3
  output->write_event( *event );
#else
  output->write_event( event.get() );
#endif
  event->clear();
}

void
HepMCHandler::fillEvent( const Event* evt )
{
  event->clear();

  // general information
#ifdef HEPMC_VERSION3
  HepMC::GenCrossSectionPtr xs = HepMC::make_shared<HepMC::GenCrossSection>();
  xs->set_cross_section( cross_sect_, cross_sect_err_ );
  event->add_attribute( "AlphaQCD", HepMC::make_shared<HepMC::DoubleAttribute>( Constants::alphaQCD ) );
  event->add_attribute( "AlphaEM", HepMC::make_shared<HepMC::DoubleAttribute>( Constants::alphaEM ) );
#else
  HepMC::GenCrossSection xs;
  xs.set_cross_section( cross_sect_, cross_sect_err_ );
  event->set_alphaQCD( Constants::alphaQCD );
  event->set_alphaQED( Constants::alphaEM );
#endif
  event->set_cross_section( xs );

  event->set_event_number( event_num_ );
  event->weights().push_back( 1. ); //FIXME we generate unweighted events

  // filling the particles content
  const HepMC::FourVector origin( 0., 0., 0., 0. );
  Particles part_vec = evt->particles();

  int cm_id = 0, idx = 1;

#ifdef HEPMC_VERSION3
  HepMC::GenVertexPtr v1 = HepMC::make_shared<HepMC::GenVertex>( origin ),
                      v2 = HepMC::make_shared<HepMC::GenVertex>( origin ),
                      vcm = HepMC::make_shared<HepMC::GenVertex>( origin );
#else
  HepMC::GenVertex* v1 = new HepMC::GenVertex( origin ),
                   *v2 = new HepMC::GenVertex( origin ),
                   *vcm = new HepMC::GenVertex( origin );
#endif

  for ( unsigned int i = 0; i < part_vec.size(); ++i ) {

    const Particle part_orig = part_vec.at( i );
    HepMC::FourVector pmom( part_orig.momentum().px(),
                            part_orig.momentum().py(),
                            part_orig.momentum().pz(),
                            part_orig.energy() );
#ifdef HEPMC_VERSION3
    HepMC::GenParticlePtr part = HepMC::make_shared<HepMC::GenParticle>( pmom, part_orig.integerPdgId(), part_orig.status() );
#else
    HepMC::GenParticle* part = new HepMC::GenParticle( pmom, part_orig.integerPdgId(), part_orig.status() );
    part->suggest_barcode( idx++ );
#endif

    const ParticlesIds moth = part_orig.mothers();

    switch ( part_orig.role() ) {
      case Particle::IncomingBeam1: { v1->add_particle_in( part ); } break;
      case Particle::IncomingBeam2: { v2->add_particle_in( part ); } break;
      case Particle::OutgoingBeam1: { v1->add_particle_out( part ); } break;
      case Particle::OutgoingBeam2: { v2->add_particle_out( part ); } break;
      case Particle::Parton1:       { v1->add_particle_out( part ); vcm->add_particle_in( part ); } break;
      case Particle::Parton2:       { v2->add_particle_out( part ); vcm->add_particle_in( part ); } break;
      case Particle::Parton3:       { v2->add_particle_out( part ); vcm->add_particle_in( part ); } break;
      case Particle::Intermediate:  { cm_id = i; continue; } break;
      case Particle::CentralSystem:
      default: {
        if ( moth.size() == 0 ) { continue; }
        if ( *moth.begin() == cm_id ) { vcm->add_particle_out( part ); }
        else {
          std::cout << "other particle!!" << std::endl;
          continue;
          //FIXME secondary products... to be implemented!
        }
      } break;
    }
    idx++;
  }
  event->add_vertex( v1 );
  event->add_vertex( v2 );
  event->add_vertex( vcm );

#ifndef HEPMC_VERSION3
  event->set_beam_particles( *v1->particles_in_const_begin(), *v2->particles_in_const_begin() );
  event->set_signal_process_vertex( *v1->vertices_begin() );
  event->set_beam_particles( *v1->particles_in_const_begin(), *v2->particles_in_const_end() );
#endif

  event_num_++;
}

#endif
