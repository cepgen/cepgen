#ifndef Jetset7Hadroniser_h
#define Jetset7Hadroniser_h

#include <algorithm>

#include "GenericHadroniser.h"

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
class Jetset7Hadroniser : public GenericHadroniser
{
 public:
  Jetset7Hadroniser();
  ~Jetset7Hadroniser();
  bool Hadronise(Particle* part_);
  bool Hadronise(Event* ev_);
 private:
  /**
   * Gives the mass for a parton/particle
   * @param[in] pdgid_ PDG id of the parton/particle
   */
  inline static double ulmass(int pdgid_) { return (double)ulmass_(pdgid_); };
  /**
   * Administrate the fragmentation and decay chain. It may be called several times, but only entries which have not yet been treated (more precisely, have 1 \f$\leq\f$ KS \f$\leq\f$ 10) can be affected by further calls. This may apply if more jets/particles have been added by the user, or if particles previously considered stable are now allowed to decay. The actions that will be taken during a LUEXEC call can be tailored extensively via the LUDAT1 - LUDAT3 commonblocks, in particular by setting the MSTJ values suitably.
   */
  inline static void luexec() { luexec_(); };
  /**
   * Set the value of any variable residing in the commmonblocks LUJETS, LUDAT1, LUDAT2, LUDAT3, LUDAT4, or LUDATR. This is done in a more controlled fashion than by directly including the commonblocks in the user program, in that array bounds are checked and the old and new values for the variable changed is written to the output for reference.
   * @param[in] line_ The line to be parsed and fed to the Jetset instance
   */
  inline static void lugive(const std::string &line_) { lugive_(line_.c_str(),line_.length()); };
  /**
   * Give the charge for a parton/particle
   * @param[in] pdgid_ PDG id of the parton/particle
   */
  inline static float luchge(int pdgid_) { return luchge_(pdgid_)/3.; };
  /**
   * The @a mlist_ parameter can take these values :
   *  - 0 : writes a header with program version number and last date of change; is mostly for internal use.
   *  - 1 : gives a simple list of current event record, in an 80 column format suitable for viewing directly on the computer terminal. For each entry, the following information is given: 
   *    - the entry number I,
   *    - the parton/particle name (see below),
   *    - the status code KS (K(I,1)),
   *    - the flavour code KF (K(I,2)),
   *    - the line number of the mother (K(I,3)), and
   *    - the three-momentum, energy and mass (P(I,1) - P(I,5)).
   *
   *    If MSTU(3) is nonzero, lines immediately after the event record proper are also listed. A final line contains information on total charge, momentum, energy and invariant mass.
   *    The particle name is given by a call to the routine LUNAME.
   *    For an entry which has decayed/fragmented (KS = 11 - 20), this particle name is given within parantheses. Similarly, a documentation line (KS = 21 - 30) has the name enclosed in expression signs (!...!) and an event/jet axis information line the name within inequality signs (<...>). If the last character of the name is a ?, it is a signal that the complete name has been truncated to fit in, and can therefore not be trusted; this is very rare.
   *    For partons which have been arranged along strings (KS = 1, 2, 11 or 12), the end of the parton name column contains information about the colour string arrangement:
   *    - a A for the first entry of a string,
   *    - an I for all intermediate ones, and
   *    - a V for the final one (a poor man's vertical rendering of the doublesided arrow <---->).
   *
   *    It is possible to insert lines just consisting of sequences of ====== to separate different sections of the event record, see MSTU(70) - MSTU(80).
   *  - 2 : gives a more extensive list of the current event record, in a 132 column format, suitable for printers or workstations.
   *    For each entry, the following information is given:
   *    - the entry number I,
   *    - the parton/particle name (with padding as described for @a mlist_ = 1),
   *    - the status code KS (K(I,1)),
   *    - the flavour code KF (K(I,2)),
   *    - the line number of the mother (K(I,3)),
   *    - the decay product/colour flow pointers (K(I,4), K(I,5)), and
   *    - the three-momentum, energy and mass (P(I,1) - P(I,5)).
   *
   *    If MSTU(3) is nonzero, lines immediately after the event record proper are also listed.
   *    A final line contains information on total charge, momentum, energy and invariant mass. Lines with only ====== may be inserted as for MLIST(1).
   *  - 3 : gives the same basic listing as = 2, but with an additional line for each entry containing information on production vertex position and time (V(I,1) - V(I,4)) and, for unstable particles, invariant lifetime (V(I,5)).
   *  - 11 : provides a simple list of all parton/particle codes defined in the program, with KF code and corresponding particle name.
   *    The list is grouped by particle kind, and only within each group in ascending order.
   *  - 12 : provides a list of all parton/particle and decay data used in the program.
   *    Each parton/particle code is represented by one line containing:
   *    - KF flavour code,
   *    - KC compressed code,
   *    - particle name,
   *    - antiparticle name (where appropriate),
   *    - electrical and colour charge (stored in KCHG),
   *    - mass,
   *    - resonance width and maximum broadening,
   *    - average invariant lifetime (in PMAS) and whether the particle is considered stable or not (in MDCY).
   *
   *    Immediately after a particle, each decay channel gets one line, containing:
   *    - decay channel number (IDC read from MDCY),
   *    - on/off switch for the channel,
   *    - matrix element type (MDME),
   *    - branching ratio (BRAT), and
   *    - decay products (KFDP).
   *
   *    The MSTU(14) flag can be used to set the maximum flavour for which particles are listed, with the default (= 0) corresponding to separately defined ones (KC > 100 if KF > 0).
   *    In order to keep the size down, decay modes of heavy hadrons collectively defined are never listed; these have KC codes 84 - 88, where the relevant information may be found.
   *  - 13 : gives a list of current parameter values for MSTU, PARU, MSTJ and PARJ, and the first 200 entries of PARF.
   *    This is useful to keep check of which default values were changed in a given run.
   * @brief List an event, jet or particle data, or current parameter values
   * @param[in] mlist_ Determines what is to be listed (see detailed description)
   */
  inline static void lulist(int mlist_) { lulist_(mlist_); };
  /**
   * Give the parton/particle name (as a character string).
   * @param[in] pdgid_ PDG id of the parton/particle
   */
  inline static std::string luname(int pdgid_) {
    char out[NAME_CHR];
    std::string s;
    luname_(pdgid_, out, NAME_CHR);
    s = std::string(out, NAME_CHR);
    s.erase(remove(s.begin(), s.end(), ' '), s.end());
    return s;
  };
  /**
   * Connect a number of previously defined partons into a string configuration.
   *
   * Initially the partons must be given with status codes (KS = K(I,1)) 1, 2 or 3.
   *
   * Afterwards the partons all have status code 3, i.e. are given with full colour flow information.
   *
   * Compared to the normal way of defining a parton system, the partons need therefore not appear in the same sequence in the event record as they are assumed to do along the string.
   *
   * It is also possible to call LUSHOW for all or some of the entries making up the string formed by @a lujoin.
   * @param[in] njoin_ Number of particles to be joined by one string
   * @param[in] ijoin_ List of particles to join in the colour flow. An one-dimensional array, of size at least @a njoin_. The @a njoin_ first numbers are the positions of the partons that are to be joined, given in the order the partons are assumed to appear along the string. If the system consists entirely of gluons, the string is closed by connecting back the last to the first entry.
   * @note Only one string (i.e. one colour singlet) may be defined per call, but one is at liberty to use any number of @a lujoin calls for a given event. The program will check that the parton configuration specified makes sense, and not take any action unless it is. Note, however, that an initially sensible parton configuration may become nonsensical, if only some of the partons are reconnected, while the others are left unchanged.
   */
  inline static void lujoin(int njoin_, int ijoin_[2]) { return lujoin_(njoin_,*ijoin_); };
  bool PrepareHadronisation(Event *ev_);
};

#endif
