#ifndef _PYTHIA6HADRONISER_H
#define _PYTHIA6HADRONISER_H

#include <algorithm>

#include "hadroniser.h"
extern "C"
{
  extern double pymass_(int&);
  extern void pyexec_();
  extern void pygive_(const char*,int);
  extern void pyckbd_();
  extern void pylist_(int&);
  extern void pyjoin_(int&,int&);
  extern void pyname_(int&,char*,int);
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
  bool Hadronise(Event* ev_);
 private:
  inline static double pymass(int pdgid_) { return pymass_(pdgid_); };
  inline static void pyexec() { pyexec_(); };
  inline static void pyckbd() { pyckbd_(); };
  inline static void pygive(const std::string &line_) { pygive_(line_.c_str(),line_.length()); };
  inline static void pylist(int mlist_) { pylist_(mlist_); };
  inline static std::string pyname(int pdgid_) {
    char out[6];
    std::string s;
    pyname_(pdgid_, out, 6);
    s = std::string(out,6);
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
  };
  /**
   * @brief Connect entries with colour flow information
   * @param njoin_ Number of particles to join in the colour flow
   * @param ijoin_ List of particles to join in the colour flow
   */
  inline static void pyjoin(int njoin_, int ijoin_[2]) { return pyjoin_(njoin_,*ijoin_); };
};

#endif
