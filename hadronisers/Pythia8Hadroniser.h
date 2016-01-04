#ifndef Pythia8Hadroniser_h
#define Pythia8Hadroniser_h

#include <algorithm>

#include "../include/GenericHadroniser.h"

/**
 * Full interface to the Pythia8 hadronisation algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
 * @brief Pythia8 hadronisation algorithm
 */
class Pythia8Hadroniser : public GenericHadroniser
{
 public:
  inline Pythia8Hadroniser() {;}
  inline ~Pythia8Hadroniser() {;}
  bool Hadronise(Particle* part_);
  bool Hadronise(Event* ev_);
 private:
  bool PrepareHadronisation(Event *ev_);
};

#endif
