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
      //GenEvent hev;
      const ::HepMC::FourVector origin(0., 0., 0., 0.);
      ::HepMC::GenVertex *v1, *v2, *DecVtx=0;

      ConstParticlesRef part = ev.GetConstParticlesRef();

      std::vector<::HepMC::GenParticle*> hepmc_particles;
      for (unsigned int j=0; j<part.size(); j++) {
        const Particle* p = part[j]; Particle::Momentum p_m = p->GetMomentum();
        ::HepMC::FourVector pmom(p_m.Px(), p_m.Py(), p_m.Pz(), p_m.E());
        ::HepMC::GenParticle* gp = new ::HepMC::GenParticle(pmom, p->GetIntPDGId(), p->status);
        gp->suggest_barcode(j+1);
        hepmc_particles.push_back(gp);

        p->Dump();

        if (j==0) { // first incoming proton
          v1 = new ::HepMC::GenVertex(origin);
          v1->add_particle_in(gp);
          last_event->add_vertex(v1);
          continue;
        }
        if (j==1) { // second incoming proton
          v2 = new ::HepMC::GenVertex(origin);
          v2->add_particle_in(gp);
          last_event->add_vertex(v2);
          continue;
        }

        ParticlesIds moth = p->GetMothersIds();
        if (moth.size()==0) {
          std::cout << "particle with no mother:" << std::endl;
          p->Dump();
          continue;
        }
	
        const int parent = *moth.begin();
        if (parent==0) { v1->add_particle_out(gp); }
        else if (parent==1) { v2->add_particle_out(gp); }
	else {
	  std::cout << "retrieving particle with barcode " << parent+1 << std::endl;
	  ::HepMC::GenParticle* parentPart = last_event->barcode_to_particle(parent+1);
	  if (!parentPart) continue; //FIXME
	  parentPart->set_status(2); // reset status, to mark that it's decayed
	  DecVtx = new ::HepMC::GenVertex(origin);
	  DecVtx->add_particle_in(parentPart);
	  DecVtx->add_particle_out(gp);
	  last_event->add_vertex(DecVtx);
	}
        for (unsigned int k=j+1; k<part.size(); k++) { // the pointer is shifted by -1, c++ style
          const Particle* p2 = part[k]; Particle::Momentum p_m2 = p2->GetMomentum();
	  std::cout << "---> retrieving next particle with id=" << p2->id << " and parent=" << *p2->GetMothersIds().begin() << std::endl; 
          if (*p2->GetMothersIds().begin()!=parent) break; // another parent particle, break the loop
          ::HepMC::FourVector pmomN(p_m2.Px(), p_m2.Py(), p_m2.Pz(), p_m2.E());
          ::HepMC::GenParticle* daughterN = new ::HepMC::GenParticle(pmomN, p2->GetIntPDGId(), p2->status); //FIXME
          daughterN->suggest_barcode(k+1);
          if (parent==1) { v1->add_particle_out(daughterN); continue; }
          if (parent==2) { v2->add_particle_out(daughterN); continue; }
          DecVtx->add_particle_out(daughterN);
	}
      }
      last_event->set_beam_particles(hepmc_particles[0], hepmc_particles[1]);
      last_event->set_signal_process_vertex(*(v1->vertices_begin()));
last_event->print();
      return last_event;
    }
  }
}
