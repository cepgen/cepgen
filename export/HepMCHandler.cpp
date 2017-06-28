#include "HepMCHandler.h"

OutputHandler::HepMCHandler::HepMCHandler( const char* filename, const ExportHandler::OutputType& type ) :
  ExportHandler( type )//, event( 0 )
{
#ifdef HEPMC_VERSION3
  output = new HepMC::WriterAscii( filename );
#else
  output = new HepMC::IO_GenEvent( filename );
#endif
}

OutputHandler::HepMCHandler::~HepMCHandler()
{
  //if ( event ) delete event;
  delete output;
}

void
OutputHandler::HepMCHandler::operator<<( const Event* evt )
{
  fillEvent( evt );
  //event->print();
#ifdef HEPMC_VERSION3
  output->write_event( *event );
#else
  *output << event;
#endif
  clearEvent();
}

void
OutputHandler::HepMCHandler::clearEvent()
{
  //if ( event ) delete event;
}

void
OutputHandler::HepMCHandler::fillEvent( const Event* evt )
{
  event = std::make_shared<HepMC::GenEvent>();
evt->dump();

  // general information
#ifdef HEPMC_VERSION3
  HepMC::GenCrossSectionPtr xs = std::make_shared<HepMC::GenCrossSection>();
  xs->set_cross_section( cross_sect_, cross_sect_err_ );
  event->set_cross_section( xs );
#else
  event->set_cross_section( cross_sect_, cross_sect_err_ );
#endif

  event->set_event_number( event_num_ );
  event->weights().push_back( 1. ); //FIXME we generate unweighted events

  // filling the particles content
  const HepMC::FourVector origin( 0., 0., 0., 0. );
  ConstParticlesRef part_vec = evt->constParticlesRef();

  int cm_id = 0, idx = 1;

#ifdef HEPMC_VERSION3
  HepMC::GenVertexPtr v1( new HepMC::GenVertex( origin ) ), v2( new HepMC::GenVertex( origin ) ), vcm( new HepMC::GenVertex( origin ) );
#else
  HepMC::GenVertex *v1 = new HepMC::GenVertex( origin ), *v2 = new HepMC::GenVertex( origin ), *vcm = new HepMC::GenVertex( origin );
#endif
  
  event->add_vertex( v1 );
  event->add_vertex( v2 );
  event->add_vertex( vcm );

  for ( unsigned int i=0; i<part_vec.size(); i++ ) {

    const Particle* part_orig = part_vec.at( i );
    HepMC::FourVector pmom( part_orig->momentum().px(),
                            part_orig->momentum().py(),
                            part_orig->momentum().pz(),
                            part_orig->energy() );
#ifdef HEPMC_VERSION3
    HepMC::GenParticlePtr part( new HepMC::GenParticle( pmom, part_orig->integerPdgId(), part_orig->status ) );
#else
    HepMC::GenParticle* part = new HepMC::GenParticle( pmom, part_orig->integerPdgId(), part_orig->status );
    part->suggest_barcode( idx++ );
#endif

    const ParticlesIds moth = part_orig->mothersIds();

    switch ( part_orig->role ) {
      case Particle::IncomingBeam1: { v1->add_particle_in( part ); } break;
      case Particle::IncomingBeam2: { v2->add_particle_in( part ); } break;
      case Particle::OutgoingBeam1: { v1->add_particle_out( part ); } break;
      case Particle::OutgoingBeam2: { v2->add_particle_out( part ); } break;
      case Particle::Parton1:       { v1->add_particle_out( part ); vcm->add_particle_in( part ); } break;
      case Particle::Parton2:       { v2->add_particle_out( part ); vcm->add_particle_in( part ); } break;
      case Particle::Parton3:       { v2->add_particle_out( part ); vcm->add_particle_in( part ); } break;
      case Particle::CentralSystem: { cm_id = i; continue; } break;
      case Particle::CentralParticle1:
      case Particle::CentralParticle2:
      default: {
        if ( moth.size()==0 ) { continue; }
        if ( *moth.begin()==cm_id ) { vcm->add_particle_out( part ); }
        else {
          std::cout << "other particle!!" << std::endl;
          continue;
          //FIXME secondary products... to be implemented!
        }
      } break;
    }
    idx++;
  }
#ifndef HEPMC_VERSION3
  event->set_beam_particles( *v1->particles_in_const_begin(), *v2->particles_in_const_begin() );
  event->set_signal_process_vertex( *v1->vertices_begin() );
  event->set_beam_particles( *v1->particles_in()_const_begin(), *v2->particles_in_const_end() );
#endif

  event_num_++;
}
