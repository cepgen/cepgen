#ifndef utils_h
#define utils_h

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <ctime>
#include <stdexcept>
#include <stdarg.h>  // For va_start, etc.

#include "Exception.h"

static std::stringstream error;

/// Provide a random number generated along a uniform distribution between 0 and 1
inline double drand() { srand (time(NULL)); return static_cast<double>(rand())/RAND_MAX; }

/// Format a string using a printf style format descriptor.
std::string Form(const std::string fmt, ...);

/**
 * A generic timer to extract the processing time between two steps in this software's flow
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
class Timer
{
 public:
  inline Timer() { clock_gettime(CLOCK_REALTIME, &beg_); }
  /**
   * Get the time elapsed since the last @a reset call (or class construction)
   * @return Elapsed time (since the last reset), in seconds
   */
  inline double elapsed() {
    clock_gettime(CLOCK_REALTIME, &end_);
    return end_.tv_sec -beg_.tv_sec+(end_.tv_nsec - beg_.tv_nsec)/1.e9;
  }
  /// Reset the clock counter
  inline void reset() {
    clock_gettime(CLOCK_REALTIME, &beg_);
  }
 private:
  /// Timestamp marking the beginning of the counter
  timespec beg_;
  /// Timestamp marking the end of the counter
  timespec end_;
};

/// List of physical constants useful that may be used for the matrix element definition
class Constants
{
 public:
  /// Electromagnetic coupling constant \f$\alpha_{em}=\frac{e^2}{4\pi\epsilon_0\hbar c}\f$
  static double AlphaEM;
  /// \f$\frac{1}{(\hbar c)^2}~[\mathrm b^{-1}]\f$?
  static double muBarn;
  /// Good ole' pi!
  static double Pi;
  /// Conversion factor between GeV^2 and barn
  static double GeV2toBarn;
  static double sconstb;
  static double AlphaReduced;
};

/**
 * Define modified variables of integration to avoid peaks integrations (see @cite Vermaseren1983347 for details)
 * Return a set of two modified variables of integration to maintain the stability of the integrant. These two new variables are :
 * - \f$y_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\f$ the new variable
 * - \f$\mathrm dy_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\log\frac{x_{min}}{x_{max}}\f$, the new variable's differential form
 * @brief Redefine the variables of integration in order to avoid the strong peaking of the integrant.
 * @param[in] expo_ Exponant
 * @param[in] xmin_ Minimal value of the variable
 * @param[in] xmax_ Maximal value of the variable
 * @param[out] out_ The new variable definition
 * @param[out] dout_ The differential variant of the new variable definition
 * @param[in] var_name_ The variable name
 * @note This method overrides the set of `mapxx` subroutines in ILPAIR, with a slight difference according to the sign of the
 *  \f$\mathrm dy_{out}\f$ parameter :
 *  - left unchanged :
 * > `mapw2`, `mapxq`, `mapwx`, `maps2`
 *  - opposite sign :
 * > `mapt1`, `mapt2`
 */
void Map(double expo_, double xmin_, double xmax_, double* out_, double* dout_, const std::string& var_name_="");
void Mapla(double,double,int,double,double,double*,double*);

/**
 * Generate random number with Breit-Wigner distribution
 * @return Random number between emin_ and emax_ with Breit-Wigner distribution:
 *  \f$\frac{1}{(E-E_r)^2+\Gamma^2/4}\f$
 * @param er_ Maximum of distribution
 * @param gamma_ Width of distribution
 * @param emin_ Minimal value of RanBW
 * @param emax_ Maximal value of RanBW
 * @date 11 Apr 2014
 */
double BreitWigner(double er, double gamma, double emin, double emax, double e=-1.);
/// Convert a polar angle to a pseudo-rapidity
inline double ThetaToEta(double theta_) { return -log(tan(theta_/180.*Constants::Pi/2.)); }
/// Convert a pseudo-rapidity to a polar angle
inline double EtaToTheta(double eta_) { return 2.*atan(exp(-eta_))*180./Constants::Pi; }
/// Convert a pseudo-rapidity to a rapidity
double EtaToY(double eta_, double m_, double pt_);

#endif
