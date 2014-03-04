#include "../include/event.h"
#include "../include/pythia6hadroniser.h"
#include "../include/jetset7hadroniser.h"
#include "../include/pythia8hadroniser.h"

int main() {

  //Event ev;
  
  Particle test1(1, 2212);
  std::cout << "--> " << test1.M() << std::endl;
  Particle test2(1, 2203);
  std::cout << "--> " << test2.M() << std::endl;
  Particle test3(1, 2101);
  std::cout << "--> " << test3.M() << std::endl;
  Particle test4(1, 2103);
  std::cout << "--> " << test4.M() << std::endl;
  Particle test5(1, 1);
  std::cout << "--> " << test5.M() << std::endl;
  Particle test6(1, 2);
  std::cout << "--> " << test6.M() << std::endl;

  //Particle pp1(1, 2212);
  //pp1.P(1., -1., 10.);
  //pp1.status = 1;
  //ev.AddParticle(&pp1);

  //Particle p1(5, 1);
  //p1.P(1., -15., 100.);
  //p1.status = 3;
  //ev.AddParticle(&p1);

  /*Particle p12(5, 2);
  p1.P(2., -25., 200.);
  p1.status = 3;
  ev.AddParticle(&p12);*/

  //Particle p2(5, 2203);
  //p2.P(10., 5., 3200.);
  //p2.status = 3;
  //ev.AddParticle(&p2);

  //ev.Dump();
  //p.E(1.312);

  //Pythia8Hadroniser py;
  //py.Hadronise(&ev);

  //Jetset7Hadroniser js;
  //js.Hadronise(&ev);
  //ev.Dump(true);

  return 0;
}
