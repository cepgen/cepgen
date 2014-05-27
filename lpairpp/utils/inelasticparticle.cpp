//#include "../include/pythia6hadroniser.h"
#include "../include/jetset7hadroniser.h"
//#include "../include/pythia8hadroniser.h"
#include "../include/physics.h"

int main() {

  //Event ev;
  
  //Particle p1(5, 1);
  //p1.P(1., -15., 100.);
  //p1.status = 3;
  //ev.AddParticle(&p1);

  //Particle p2(5, 2203);
  //p2.P(10., 5., 3200.);
  //p2.status = 3;
  //ev.AddParticle(&p2);

  //Jetset7Hadroniser js;
  //js.Hadronise(&ev);
  //ev.Dump(true);

  Particle p1(1, 11);
  Particle p2(2, 2212);
  p1.P(0., 0., 27.55);
  p2.P(0., 0., -820.);

  PhysicsBoundaries pb;

  Particles epa = EPA(p1, p2, 1, pb);

  Particles::iterator part;
  for (part=epa.begin(); part!=epa.end(); part++) {
    part->Dump();
  }

  return 0;
}
