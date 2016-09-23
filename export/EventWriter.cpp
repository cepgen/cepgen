#include "EventWriter.h"

EventWriter::EventWriter(const OutputType& type, const char* filename) :
  fType(type),
  fHepMCOutput(0)
{
  switch (fType) {
    case HepMC: {
      fHepMCOutput = new HepMC::IO_GenEvent(filename);
    } break;
    default: return;
  }
}

EventWriter::~EventWriter()
{
  // HepMC persistent objects
  //if (fHepMCEvent) delete fHepMCEvent;
  if (fHepMCOutput) delete fHepMCOutput;
}

void
EventWriter::operator<<(const Event* evt)
{
  switch (fType) {
    case HepMC: {
      const HepMC::GenEvent* ev = getHepMCEvent(evt);
      (*fHepMCOutput) << ev;
      delete ev;
    } break;
    default: return;
  }
}

HepMC::GenParticle*
EventWriter::getHepMCParticle(const Particle* part) const
{
  HepMC::FourVector pmom(part->GetMomentum().Px(),
                         part->GetMomentum().Py(),
                         part->GetMomentum().Pz(),
                         part->E());
  return new HepMC::GenParticle(pmom, part->GetIntPDGId(), part->status);
}

HepMC::GenEvent*
EventWriter::getHepMCEvent(const Event* evt) const
{
  HepMC::GenEvent* out = new HepMC::GenEvent;
  const HepMC::FourVector origin(0., 0., 0., 0.);
  ConstParticlesRef part_vec = evt->GetConstParticlesRef();

  int cm_id = 0, idx = 1;

  HepMC::GenVertex *v1 = new HepMC::GenVertex(origin), 
                   *v2 = new HepMC::GenVertex(origin),
                   *vcm = new HepMC::GenVertex(origin);

  out->add_vertex(v1);
  out->add_vertex(v2);
  out->add_vertex(vcm);

  for (unsigned int j=0; j<part_vec.size(); j++) {

    HepMC::GenParticle* part = getHepMCParticle(part_vec[j]);
    const ParticlesIds moth = part_vec[j]->GetMothersIds();

    switch (part_vec[j]->role) {
      case Particle::IncomingBeam1: {
        part->suggest_barcode(idx++);
        v1->add_particle_in(part);
      } break;
      case Particle::IncomingBeam2: {
        part->suggest_barcode(idx++);
        v2->add_particle_in(part);
      } break;
      case Particle::OutgoingBeam1: {
        part->suggest_barcode(idx++);
        v1->add_particle_out(part);
      } break;
      case Particle::OutgoingBeam2: {
        part->suggest_barcode(idx++);
        v2->add_particle_out(part);
      } break;
      case Particle::Parton1: {
        part->suggest_barcode(idx++);
        v1->add_particle_out(part);
        vcm->add_particle_in(part);
      } break;
      case Particle::Parton2: {
        part->suggest_barcode(idx++);
        v2->add_particle_out(part);
        vcm->add_particle_in(part);
      } break;
      case Particle::Parton3: {
        part->suggest_barcode(idx++);
        v2->add_particle_out(part);
        vcm->add_particle_in(part);
      } break;
      case Particle::CentralSystem: {
        cm_id = j; continue;
      } break;
      case Particle::CentralParticle1:
      case Particle::CentralParticle2:
      default: {
        part->suggest_barcode(idx++);
        if (moth.size()==0) { continue; }
        if (*moth.begin()==cm_id) { vcm->add_particle_out(part); }
        else {
          std::cout << "other particle!!" << std::endl;
          continue;
          //FIXME secondary products... to be implemented!
        }
      } break;
    }
  }


  out->set_beam_particles(*v1->particles_in_const_begin(), *v2->particles_in_const_begin());
  out->set_signal_process_vertex(*(v1->vertices_begin()));

  HepMC::GenCrossSection xs;
  xs.set_cross_section(fCrossSect, fCrossSectErr);
  out->set_cross_section(xs);

  out->set_event_number(fEventNum);

  return out;
}
