#include "../include/inelastic.h"

int main() {
  
  InelasticParticle p;

  p.pdgId = 2;
  p.M(1.312);
  //std::cout << p.P(10., 5., 3200.) << std::endl;
  //std::cout << p.M() << std::endl;
  //p.Dump();
  std::cout << "hadronisation : " << p.Hadronise("pythia6") << std::endl;

  return 0;
}
