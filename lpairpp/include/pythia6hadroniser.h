#ifndef _PYTHIA6HADRONISER_H
#define _PYTHIA6HADRONISER_H

#include <algorithm>

#include "hadroniser.h"

#define NAME_CHR 16

extern "C"
{
  extern double pymass_(int&);
  extern void pyexec_();
  extern void pygive_(const char*,int);
  extern void pyckbd_();
  extern void pylist_(int&);
  extern void pyjoin_(int&,int&);
  extern void pyname_(int&,char*,int);
  extern double pyp_(int&,int&);
  extern struct
  {
    int n, npad, k[5][4000];
    double p[5][4000], v[5][4000];
  } pyjets_;
}

/**
 * Full interface to the Pythia6 @cite Sjostrand:2006za algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
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
  inline static double pyp(int role_, int qty_) { return pyp_(role_,qty_); };
  inline static std::string pyname(int pdgid_) {
    char out[NAME_CHR];
    std::string s;
    pyname_(pdgid_, out, NAME_CHR);
    s = std::string(out, NAME_CHR);
    //s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    s.erase(remove(s.begin(), s.end(), ' '), s.end());
    return s;
  };
  /**
   * @brief Connect entries with colour flow information
   * @param[in] njoin_ Number of particles to join in the colour flow
   * @param[in] ijoin_ List of particles unique identifier to join in the colour flow
   */
  inline static void pyjoin(int njoin_, int ijoin_[2]) { return pyjoin_(njoin_,*ijoin_); };
  void PrepareHadronisation(Event *ev_);
};

#endif
