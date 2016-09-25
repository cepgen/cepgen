#include "HepMCHandler.h"

OutputHandler::HepMCHandler::HepMCHandler(const char* filename):
  event(0)
{
#ifndef HEPMC_VERSION_CODE
  output = new HepMC::IO_GenEvent(filename);
#else
  output = new HepMC::WriterAscii(filename);
#endif
}

OutputHandler::HepMCHandler::~HepMCHandler()
{
  delete output;
}

void
OutputHandler::HepMCHandler::operator<<(const Event* evt)
{
  fillEvent(evt);
#ifndef HEPMC_VERSION_CODE
  *output << event;
#else
  output->write_event(*event);
#endif
  clearEvent();
}

void
OutputHandler::HepMCHandler::clearEvent()
{
  if (event) delete event;
}

void
OutputHandler::HepMCHandler::fillEvent(const Event* evt)
{
  event = new HepMC::GenEvent;

  // general information
  HepMC::GenCrossSection xs;
  xs.set_cross_section(fCrossSect, fCrossSectErr);
#ifndef HEPMC_VERSION_CODE
  event->set_cross_section(xs);
#else
  event->set_cross_section(HepMC::GenCrossSectionPtr(&xs));
#endif

  event->set_event_number(fEventNum);

  // filling the particles content
  const HepMC::FourVector origin(0., 0., 0., 0.);
  ConstParticlesRef part_vec = evt->GetConstParticlesRef();

  int cm_id = 0, idx = 1;

#ifndef HEPMC_VERSION_CODE
  HepMC::GenVertex *v1 = new HepMC::GenVertex(origin),
                   *v2 = new HepMC::GenVertex(origin),
                   *vcm = new HepMC::GenVertex(origin);
#else
  HepMC::GenVertexPtr v1(new HepMC::GenVertex(origin)),
                      v2(new HepMC::GenVertex(origin)),
                      vcm(new HepMC::GenVertex(origin));
#endif
  
  event->add_vertex(v1);
  event->add_vertex(v2);
  event->add_vertex(vcm);

  for (unsigned int i=0; i<part_vec.size(); i++) {

    const Particle* part_orig = part_vec.at(i);
    HepMC::FourVector pmom(part_orig->GetMomentum().Px(),
                           part_orig->GetMomentum().Py(),
                           part_orig->GetMomentum().Pz(),
                           part_orig->E());
#ifndef HEPMC_VERSION_CODE
    HepMC::GenParticle* part = new HepMC::GenParticle(pmom, part_orig->GetIntPDGId(), part_orig->status);
    part->suggest_barcode(idx++);
#else
    HepMC::GenParticlePtr part(new HepMC::GenParticle(pmom, part_orig->GetIntPDGId(), part_orig->status));
#endif

    const ParticlesIds moth = part_orig->GetMothersIds();

    switch (part_orig->role) {
      case Particle::IncomingBeam1: { v1->add_particle_in(part); } break;
      case Particle::IncomingBeam2: { v2->add_particle_in(part); } break;
      case Particle::OutgoingBeam1: { v1->add_particle_out(part); } break;
      case Particle::OutgoingBeam2: { v2->add_particle_out(part); } break;
      case Particle::Parton1:       { v1->add_particle_out(part); vcm->add_particle_in(part); } break;
      case Particle::Parton2:       { v2->add_particle_out(part); vcm->add_particle_in(part); } break;
      case Particle::Parton3:       { v2->add_particle_out(part); vcm->add_particle_in(part); } break;
      case Particle::CentralSystem: { cm_id = i; continue; } break;
      case Particle::CentralParticle1:
      case Particle::CentralParticle2:
      default: {
        if (moth.size()==0) { continue; }
        if (*moth.begin()==cm_id) { vcm->add_particle_out(part); }
        else {
          std::cout << "other particle!!" << std::endl;
          continue;
          //FIXME secondary products... to be implemented!
        }
      } break;
    }
    idx++;
  }
#ifndef HEPMC_VERSION_CODE
  event->set_beam_particles(*v1->particles_in_const_begin(), *v2->particles_in_const_begin());
  event->set_signal_process_vertex(*v1->vertices_begin());
  event->set_beam_particles(*v1->particles_in().begin(), *v2->particles_in().end());
#endif
}
