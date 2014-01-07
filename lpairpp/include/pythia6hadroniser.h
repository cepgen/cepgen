#include "hadroniser.h"

extern "C"
{
  extern double pymass_(int&);
  extern void pyexec_();
  extern struct
  {
    int n, npad, k[5][4000];
    double p[5][4000], v[5][4000];
  } pyjets_;
}

/**
 * @brief Pythia6 hadronisation algorithm
 */
class Pythia6Hadroniser : public Hadroniser
{
 public:
  Pythia6Hadroniser();
  ~Pythia6Hadroniser();
  bool Hadronise(Particle* part_);
 private:
  inline static double pymass(int pdgid) { return pymass_(pdgid); };
  inline static void pyexec() { return pyexec_(); };
};
