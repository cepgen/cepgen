#include "EventWriter.h"

namespace EventWriter
{
  namespace HepMC
  {
    ::HepMC::GenEvent*
    Event(const ::Event& ev)
    {
      if (last_event) { delete last_event; }
      last_event = new ::HepMC::GenEvent;
      const ::HepMC::FourVector origin(0., 0., 0., 0.);
      ::HepMC::GenVertex *v1, *v2, *vcm;

      ConstParticlesRef part = ev.GetConstParticlesRef();

      std::vector<::HepMC::GenParticle*> hepmc_particles;
      int cm_id = 0, idx = 1;
      for (unsigned int j=0; j<part.size(); j++) {
        const Particle* p = part[j]; Particle::Momentum p_m = p->GetMomentum();
        ::HepMC::FourVector pmom(p_m.Px(), p_m.Py(), p_m.Pz(), p_m.E());
        ::HepMC::GenParticle* gp = new ::HepMC::GenParticle(pmom, p->GetIntPDGId(), p->status);
        hepmc_particles.push_back(gp);

        switch (p->role) {
          case Particle::IncomingBeam1: {
            gp->suggest_barcode(idx++);
            v1 = new ::HepMC::GenVertex(origin);
            v1->add_particle_in(gp);
            last_event->add_vertex(v1);
          } break;
          case Particle::IncomingBeam2: {
            gp->suggest_barcode(idx++);
            v2 = new ::HepMC::GenVertex(origin);
            v2->add_particle_in(gp);
            last_event->add_vertex(v2);
          } break;
          case Particle::OutgoingBeam1: {
            gp->suggest_barcode(idx++);
            v1->add_particle_out(gp);
          } break;
          case Particle::OutgoingBeam2: {
            gp->suggest_barcode(idx++);
            v2->add_particle_out(gp);
          } break;
          case Particle::Parton1: {
            gp->suggest_barcode(idx++);
            v1->add_particle_out(gp);
            vcm = new ::HepMC::GenVertex(origin);
            vcm->add_particle_in(gp);
            last_event->add_vertex(vcm);
          } break;
          case Particle::Parton2: {
            gp->suggest_barcode(idx++);
            v2->add_particle_out(gp);
            vcm->add_particle_in(gp);
          } break;
          case Particle::Parton3: {
            gp->suggest_barcode(idx++);
            v2->add_particle_out(gp);
            vcm->add_particle_in(gp);
          } break;
          case Particle::CentralSystem: {
            cm_id = j;
          } break;
          case Particle::CentralParticle1:
          case Particle::CentralParticle2:
          default: {
            gp->suggest_barcode(idx++);
            const ParticlesIds moth = p->GetMothersIds();
            if (moth.size()==0) { continue; }
            if (*moth.begin()==cm_id) {
              vcm->add_particle_out(gp);
            }
            else {
              //FIXME secondary products... to be implemented!
            }
          } break;
        }
      }
      last_event->set_beam_particles(hepmc_particles[0], hepmc_particles[1]);
      last_event->set_signal_process_vertex(*(v1->vertices_begin()));
      return last_event;
    }
  }
}
