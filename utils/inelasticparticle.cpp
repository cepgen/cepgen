#include "../hadronisers/Pythia6Hadroniser.h"
#include "../hadronisers/Jetset7Hadroniser.h"
//#include "../hadronisers/Pythia8Hadroniser.h"
#include "../include/Physics.h"

int main() {

  /*Event ev;
  
  Particle p1(5, 1);
  p1.P(1., -15., 100.);
  p1.status = 3;
  ev.AddParticle(&p1);

  Particle p2(5, 2203);
  p2.P(10., 5., 3200.);
  p2.status = 3;
  ev.AddParticle(&p2);*/

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

  Logger::GetInstance()->Level = Logger::Debug;
  Event ev;
  Particle q(1, Particle::uQuark);
  q.status = -2;
  q.P(0.5, 0.8, 3000.);
  ev.AddParticle(q);
  ev.Dump();
  Pythia6Hadroniser h;
  //Jetset7Hadroniser h;
  try { h.Hadronise(&ev); } catch (Exception &e) { e.Dump(); }
  ev.Dump();
  
  return 0;
}
