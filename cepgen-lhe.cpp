#include <iostream>

#include "include/MCGen.h"
#include "export/EventWriter.h"

using namespace std;

/**
 * Main caller for this Monte Carlo generator. Loads the configuration files'
 * variables if set as an argument to this program, else loads a default
 * "LHC-like" configuration, then launches the cross-section computation and
 * the events generation.
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main(int argc, char* argv[]) {
  MCGen mg;
  
  if (argc==1) InError("No config file provided.");

  Debugging(Form("Reading config file stored in %s", argv[1]));
  if (!mg.parameters->ReadConfigFile(argv[1])) {
    Information(Form("Error reading the configuration!\n\t"
                     "Please check your input file (%s)", argv[1]));
    return -1;
  }

  // We might want to cross-check visually the validity of our run
  mg.parameters->Dump();

  // Let there be cross-section...
  double xsec, err;
  mg.ComputeXsection(&xsec, &err);

  HepMC::GenCrossSection xs;
  xs.set_cross_section(xsec, err);

  if (!mg.parameters->generation) return 0;

  EventWriter::HepMC::output output("example.dat", std::ios::out);

  // The events generation starts here !
  for (int i=0; i<mg.parameters->maxgen; i++) {
    if (i%10000==0)
      cout << "Generating event #" << i+1 << endl;
    const Event ev = *mg.GenerateOneEvent();
    ev.Dump();
    HepMC::GenEvent* hev = EventWriter::HepMC::Event(ev);
    hev->set_cross_section(xs);                                                                                                                                                                                                                                      
    hev->set_event_number(i); 

    /*std::vector<HepMC::GenParticle*> hepmc_particles;
    for (int iprt=0; iprt<2; iprt++) {
      Particle* p = part[iprt]; Particle::Momentum p_m = p->GetMomentum();
      HepMC::FourVector pmom(p_m.Px(), p_m.Py(), p_m.Pz(), p_m.E());
      HepMC::GenParticle* primary = new HepMC::GenParticle(pmom, p->GetIntPDGId(), p->status); //FIXME
      primary->suggest_barcode(iprt+1);
      hepmc_particles.push_back(primary);
    }
    HepMC::GenVertex* v1 = new HepMC::GenVertex(HepMC::FourVector(0., 0., 0., 0.));
    HepMC::GenVertex* v2 = new HepMC::GenVertex(HepMC::FourVector(0., 0., 0., 0.));
    v1->add_particle_in(hepmc_particles[0]);
    v2->add_particle_in(hepmc_particles[1]);
    hev->add_vertex(v1);
    hev->add_vertex(v2);
    hev->set_beam_particles(hepmc_particles[0], hepmc_particles[1]);
    //
    HepMC::GenVertex* DecVtx = 0;
    for (int iprt=2; iprt<part.size(); iprt++) {
      Particle* p = part[iprt]; Particle::Momentum p_m = p->GetMomentum();
      int parent = *p->GetMothersIds().begin();
      HepMC::FourVector pmom(p_m.Px(), p_m.Py(), p_m.Pz(), p_m.E());
      HepMC::GenParticle* daughter = new HepMC::GenParticle(pmom, p->GetIntPDGId(), p->status); //FIXME
      daughter->suggest_barcode(iprt+1);
      if (parent==0) {
        HepMC::GenVertex* v21 = new HepMC::GenVertex(HepMC::FourVector(0., 0., 0., 0.));
        v21->add_particle_in(daughter);
        hev->add_vertex(v21);
      }
      else if (parent==1) { v1->add_particle_out(daughter); }
      else if (parent==2) { v2->add_particle_out(daughter); }
      else {
        HepMC::GenParticle* parentPart = hev->barcode_to_particle(parent);
        parentPart->set_status(2); // reset status, to mark that it's decayed
        DecVtx = new HepMC::GenVertex(HepMC::FourVector(0., 0., 0., 0.));
        DecVtx->add_particle_in(parentPart);
        DecVtx->add_particle_out(daughter);
        hev->add_vertex(DecVtx);
      }
      int iprt1;
      for (iprt1=iprt+1; iprt1<part.size(); iprt1++) { // the pointer is shifted by -1, c++ style
        Particle* p = part[iprt]; Particle::Momentum p_m = p->GetMomentum();
        if (*p->GetMothersIds().begin()!=parent) break; // another parent particle, break the loop
        HepMC::FourVector pmomN(p_m.Px(), p_m.Py(), p_m.Pz(), p_m.E());
        HepMC::GenParticle* daughterN = new HepMC::GenParticle(pmomN, p->GetIntPDGId(), p->status); //FIXME
        daughterN->suggest_barcode(iprt1+1);
        if (parent==1) { v1->add_particle_out(daughterN); }
        else if (parent==2) { v2->add_particle_out(daughterN); }
        else { DecVtx->add_particle_out(daughterN); }
      }
      iprt = iprt1-1; // reset counter such that it doesn't go over the same child more than once
      // don't forget to offset back into c++ counting, as it's already +1 forward
    }
    hev->set_signal_process_vertex(*(v1->vertices_begin()));*/


    //ev.Dump();
    output << (hev);
  }

  //mg.parameters->StoreConfigFile("lastrun.card");

  return 0;
}

