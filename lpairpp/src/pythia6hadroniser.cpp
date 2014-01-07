#include "pythia6hadroniser.h"

Pythia6Hadroniser::Pythia6Hadroniser()
{
  _name = "Pythia6";
}

Pythia6Hadroniser::~Pythia6Hadroniser()
{
}

bool
Pythia6Hadroniser::Hadronise(Particle *part_)
{
  pyjets_.p[0][0] = part_->px;
  pyjets_.p[1][0] = part_->py;
  pyjets_.p[2][0] = part_->pz;
  pyjets_.p[3][0] = part_->E();
  pyjets_.p[4][0] = part_->M();

  pyjets_.k[0][0] = 1; // status
  pyjets_.k[1][0] = 2; // particle id
  pyjets_.k[2][0] = 5; // mother
  pyjets_.k[3][0] = 0; // daughter 1
  pyjets_.k[4][0] = 0; // daughter 2

  this->pyexec();
  //pyjets_.v[0][0] = 0;
  std::cout << "[Pythia6Hadroniser::Hadronise] INFO" << std::endl;
  //part_->Dump();
  return true;
}
