#ifndef _JETSET7HADRONISER_H
#define _JETSET7HADRONISER_H

#include "hadroniser.h"

extern "C"
{
  extern double ulmass_(int&);
  extern void luexec_();
  extern void lulist_(int&);
  extern void lujoin_(int&,int&);
  extern void luname_(int&,char*,int);
  extern struct
  {
    int n, k[5][4000];
    double p[5][4000], v[5][4000];
  } lujets_;
}

/**
 * @brief Jetset7 hadronisation algorithm
 */
class Jetset7Hadroniser : public Hadroniser
{
 public:
  Jetset7Hadroniser();
  ~Jetset7Hadroniser();
  bool Hadronise(Particle* part_);
  bool Hadronise(Event* ev_);
 private:
  inline static double ulmass(int pdgid_) { return ulmass_(pdgid_); };
  inline static void luexec() { luexec_(); };
  inline static void lulist(int mlist_) { lulist_(mlist_); };
  inline static std::string luname(int pdgid_) { char out[5]; luname_(pdgid_, out, 5); return std::string(out,5); };
  /**
   * @brief Connect entries with colour flow information
   * @param njoin_ Number of particles to join in the colour flow
   * @param ijoin_ List of particles to join in the colour flow
   */
  inline static void lujoin(int njoin_, int ijoin_[2]) { return lujoin_(njoin_,*ijoin_); };
};

#endif
