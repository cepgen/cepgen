#ifndef _PYTHIA8HADRONISER_H
#define _PYTHIA8HADRONISER_H

#include "hadroniser.h"

#include "../external/pythia8175/include/Pythia.h" //FIXME FIXME FIXME FIXME
//#include "Pythia.h" //FIXME FIXME FIXME FIXME

class Pythia8Hadroniser : public Hadroniser
{
 public:
  Pythia8Hadroniser();
  ~Pythia8Hadroniser();
  bool Hadronise(Event* ev_);
};


#endif
