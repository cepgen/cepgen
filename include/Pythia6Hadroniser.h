#ifndef Pythia6Hadroniser_h
#define Pythia6Hadroniser_h

#include <algorithm>

#include "GenericHadroniser.h"

/** @brief Maximal number of characters to fetch for the particle's name */
#define NAME_CHR 16

extern "C"
{
  /** @brief Get the particle's mass in GeV from the Pythia6 module */
  extern double pymass_(int&);
  /** @brief Get the resonant particle's width in GeV from the Pythia6 module */
  //extern double pywidt_(int&);
  /** @brief Launch the Pythia6 fragmentation */
  extern void pyexec_();
  /** @brief Set a parameter value to the Pythia6 module */
  extern void pygive_(const char*,int);
  extern void pyckbd_();
  /** @brief Lists all the particles in the event in a human-readable format */
  extern void pylist_(int&);
  /** @brief Joins two coloured particles in a colour singlet */
  extern void pyjoin_(int&,int&);
  /** @brief Get a particle's human-readable name from the Pythia6 module */
  extern void pyname_(int&,char*,int);
  /** @brief Get information on a particle from the Pythia6 module */
  extern double pyp_(int&,int&);
  /** @brief Stores one parton/particle in the PYJETS common block */
  extern void py1ent_(int&,int&,double&,double&,double&);

  /** @brief Particles content of the event */
  extern struct
  {
    /** @brief Number of particles in the event */
    int n;
    int npad;
    /** @brief Particles' general information (status, PDG id, mother, daughter 1, daughter 2) */
    int k[5][4000];
    /** @brief Particles' kinematics, in GeV (px, py, pz, E, M) */
    double p[5][4000];
    /** @brief Primary vertex for the particles */
    double v[5][4000];
  } pyjets_;
}

/**
 * Full interface to the Pythia6 @cite Sjostrand:2006za algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
 * @brief Pythia6 hadronisation algorithm
 */
class Pythia6Hadroniser : public GenericHadroniser
{
 public:
  Pythia6Hadroniser();
  ~Pythia6Hadroniser();
  bool Hadronise(Particle* part_);
  bool Hadronise(Event* ev_);
 private:
  inline static double pymass(int pdgid_) { return pymass_(pdgid_); };
  //inline static double pywidt(int pdgid_) { return pywidt_(pdgid_); };
  inline static void pyexec() { pyexec_(); };
  inline static void pyckbd() { pyckbd_(); };
  inline static void pygive(const std::string &line_) { pygive_(line_.c_str(),line_.length()); };
  inline static void pylist(int mlist_) { pylist_(mlist_); };
  inline static double pyp(int role_, int qty_) { return pyp_(role_,qty_); };
  //inline static void py1ent(int* kf_, double* pe_, double theta_, double phi_) { int one=1; py1ent_(&one, kf_, pe_, theta_, phi_); }
  //inline static void py1ent(int* kf_, double* pe_, double theta_, double phi_) { py1ent_(1, kf_, pe_, theta_, phi_); }
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
  bool PrepareHadronisation(Event *ev_);
};

#endif
