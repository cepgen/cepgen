#ifndef Pythia8Hadroniser_h
#define Pythia8Hadroniser_h

#include <algorithm>
#include <Pythia.h>

#include "hadronisers/GenericHadroniser.h"

/**
 * Full interface to the Pythia8 hadronisation algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
 * @brief Pythia8 hadronisation algorithm
 */
class Pythia8Hadroniser : public GenericHadroniser
{
 public:
  Pythia8Hadroniser();
  ~Pythia8Hadroniser();
  bool Hadronise(Particle* part_);
  bool Hadronise(Event* ev_);
 private:
  bool PrepareHadronisation(Event *ev_);
  Pythia8::Pythia* fPy;
};

#endif
