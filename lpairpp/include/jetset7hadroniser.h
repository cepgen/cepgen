#ifndef _JETSET7HADRONISER_H
#define _JETSET7HADRONISER_H

#include <algorithm>

#include "hadroniser.h"

#define NAME_CHR 16

extern "C"
{
  extern float ulmass_(int&);
  extern void luexec_();
  extern void lugive_(const char*,int);
  extern void lulist_(int&);
  extern void lujoin_(int&,int&);
  extern void luname_(int&,char*,int);
  extern int luchge_(int&);
  extern struct
  {
    int n, k[5][4000];
    float p[5][4000], v[5][4000];
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
  inline static double ulmass(int pdgid_) { return (double)ulmass_(pdgid_); };
  inline static void luexec() { luexec_(); };
  inline static void lugive(const std::string &line_) { lugive_(line_.c_str(),line_.length()); };
  inline static float luchge(int pdgid_) { return luchge_(pdgid_)/3.; };
  inline static void lulist(int mlist_) { lulist_(mlist_); };
  inline static std::string luname(int pdgid_) {
    char out[NAME_CHR];
    std::string s;
    luname_(pdgid_, out, NAME_CHR);
    s = std::string(out, NAME_CHR);
    s.erase(remove(s.begin(), s.end(), ' '), s.end());
    return s;
  };
  /**
   * @brief Connect entries with colour flow information
   * @param njoin_ Number of particles to join in the colour flow
   * @param ijoin_ List of particles to join in the colour flow
   */
  inline static void lujoin(int njoin_, int ijoin_[2]) { return lujoin_(njoin_,*ijoin_); };
  bool PrepareHadronisation(Event *ev_);
};

#endif
