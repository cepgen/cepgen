#include "../export/EventWriter.h"

using namespace std;

int main() {

  EventWriter writer(EventWriter::HepMC, "example.dat");
  writer.SetCrossSection(1., 2.);

  Event ev;
  
  Particle p1(Particle::IncomingBeam1, Particle::Proton);
  p1.SetMomentum(1., -15., 100.);
  p1.status = Particle::Incoming;
  ev.AddParticle(p1);

  Particle p2(Particle::IncomingBeam2, Particle::Electron);
  p2.SetMomentum(10., 5., 3200.);
  p2.status = Particle::Incoming;
  ev.AddParticle(p2);

  ev.Dump();
  
  writer << &ev;

  return 0;

  //Jetset7Hadroniser js;
  //js.Hadronise(&ev);
  //ev.Dump(true);

  /*Particle *p1, *p2;
  
  p1 = new Particle(1, ELECTRON);
  p2 = new Particle(2, PROTON);
  p1->P(0., 0., 27.55);
  p2->P(0., 0., -820.);

  PhysicsBoundaries pb;
  double q2 = 1.;

  Particles epa = EPA(p1, p2, 1, pb, &q2);

  Particles::iterator part;
  for (part=epa.begin(); part!=epa.end(); part++) {
    part->Dump();
  }

  delete p1;
  delete p2;*/
  
  /*Particle vm(1, Particle::JPsi);
  vm.P(1., 1., 0.);
  try { vm.Dump(); } catch (Exception& e) { e.Dump(); }*/

  /*Logger::GetInstance()->Level = Logger::Debug;
  Event ev;
  Particle q(1, Particle::uQuark);
  q.status = -2;
  q.P(0.5, 0.8, 3000.);
  ev.AddParticle(q);
  ev.Dump();
  Pythia6Hadroniser h;
  //Jetset7Hadroniser h;
  try { h.Hadronise(&ev); } catch (Exception &e) { e.Dump(); }
  ev.Dump();*/
  
  /*Particle ele(1, Particle::Electron);
  ele.SetMomentum(0., 0., 27.5);
  //ele.SetM();
  ele.Dump();
  
  Particle pro(2, Particle::Proton);
  pro.SetMomentum(0., 0., -920.0);
  pro.Dump();
  
  cout << CMEnergy(ele, pro) << endl;*/
  
  return 0;
}
