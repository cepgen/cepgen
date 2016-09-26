#ifndef _HERWIG6HADRONISER_H
#define _HERWIG6HADRONISER_H

#include <algorithm>

#include "core/GenericHadroniser.h"

#define NMXHEP 4000
//#define NAME_CHR 16

extern "C"
{
  void hwdhad_();
  void hwaend_();
  extern struct {
    int nevhep, nhep, isthep[NMXHEP], idhep[NMXHEP];
    int jmohep[NMXHEP][2], jdahep[NMXHEP][2];
    double phep[NMXHEP][5], vhep[NMXHEP][4];
  } hepevt_;
  
  /* COMMON/FFS/TB,BT
     COMMON/SFF/IT1,IB1,IT2,IB2 */
  /*struct {
    double tb, bt;
  } ffs_;
  struct {
    int it1, ib1, it2, ib2;
  } sff_;*/
}

/**
 * @brief Herwig6 hadronisation algorithm
 */
class Herwig6Hadroniser : public GenericHadroniser
{
 public:
  Herwig6Hadroniser();
  ~Herwig6Hadroniser();
  bool Hadronise(Event* ev_);
 private:
  inline static void hwdhad() { hwdhad_(); };
  /*inline static double ulmass(int pdgid_) { return (double)ulmass_(pdgid_); };
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
  };*/
  /**
   * @brief Connect entries with colour flow information
   * @param njoin_ Number of particles to join in the colour flow
   * @param ijoin_ List of particles to join in the colour flow
   */
  //inline static void lujoin(int njoin_, int ijoin_[2]) { return lujoin_(njoin_,*ijoin_); };
};

#endif
