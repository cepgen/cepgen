#include "../include/event.h"
#include "../include/pythia6hadroniser.h"
#include "../include/jetset7hadroniser.h"

int main() {

  Event ev;
  
  Particle pp1(1, 2212);
  pp1.P(1., -1., 10.);
  pp1.status = 1;
  ev.AddParticle(&pp1);

  Particle p1(5, 1);
  p1.P(1., -15., 100.);
  p1.status = 3;
  ev.AddParticle(&p1);

  /*Particle p12(5, 2);
  p1.P(2., -25., 200.);
  p1.status = 3;
  ev.AddParticle(&p12);*/

  Particle p2(5, 2203);
  p2.P(10., 5., 3200.);
  p2.status = 3;
  ev.AddParticle(&p2);

  ev.Dump();
  //p.E(1.312);

  Pythia6Hadroniser py;
  py.Hadronise(&ev);

  //Jetset7Hadroniser js;
  //js.Hadronise(&ev);
  ev.Dump(true);

  return 0;
}
