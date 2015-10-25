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

/**
 * @brief Provides a random number generated along a uniform distribution
 * between 0 and 1
 */
inline double drand() { srand (time(NULL)); return static_cast<double>(rand())/RAND_MAX; }

/// Formats a string using a printf style format descriptor.
std::string Form(const std::string fmt, ...);

/**
 * An object which enables to extract the processing time between two steps in
 * this software's flow
 */
class Timer
{
 public:
  inline Timer() { clock_gettime(CLOCK_REALTIME, &beg_); }
  /**
   * Get the time elapsed since the last @a reset call (or class construction)
   * @return The elapsed time in seconds
   */
  inline double elapsed() {
    clock_gettime(CLOCK_REALTIME, &end_);
    return end_.tv_sec -beg_.tv_sec+(end_.tv_nsec - beg_.tv_nsec)/1000000000.;
  }
  /**
   * @brief Resets the clock counter
   */
  inline void reset() {
    clock_gettime(CLOCK_REALTIME, &beg_);
  }
 private:
  /** @brief Timestamp marking the beginning of the counter */
  timespec beg_;
  /** @brief Timestamp marking the end of the counter */
  timespec end_;
};

/** @brief Electromagnetic coupling constant \f$\alpha_{em}=\frac{e^2}{4\pi\epsilon_0\hbar c}\f$ */
#define alphaF 1./137.04
/** @brief \f$\frac{1}{(\hbar c)^2}~[\mathrm b^{-1}]\f$? */
#define muBarn 1./389.39
#define pi 3.1415926535897932384626434
#define sconst 3.89351824E8
#define sconstb 2.1868465E10
#define alphared 1.16140981417e-3

/**
 * @brief Defines modified variables of integration to avoid peaks integrations (see @cite Vermaseren1983347 for details)
 * Returns a set of two modified variables of integration to maintain the stability of the integrant. These two new variables are :
 * - \f$y_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\f$ the new variable
 * - \f$\mathrm dy_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\log\frac{x_{min}}{x_{max}}\f$, the new variable's differential form
 * @brief Redefines the variables of integration in order to avoid the strong peaking of the integrant.
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
//void Symmetrise(double, double, double*, double*);

/**
 * Lorentz boost of a 4-vector (from CERNLIB)
 * @param pi_ Input 4-vector to boost
 * @param pf_ Output boosted 4-vector
 * @author L. Pape
 * @date 20 Aug 1975
 * @author Ian McLaren (mclareni), CERN/CN
 * @date 14 Feb 1996
 */
void Lorenb(double u_, double ps_[], double pi_[], double pf_[]);

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
double RanBW(double er_, double gamma_, double emin_, double emax_);
double GenerT(double tmin_, double tmax_, double b_, double anexp_);

/**
 * Generate \f$t\f$ between @a tmin_ and @a tmax_ with a distribution according to
 * \f$\frac{e^{-bt}}{(1+t/0.71)^n}\f$
 * This is the result of the flux factor given by DONNACHIE and LANDSHOFF.
 * @note
 *  - @a t, @a tmin_, and @a tmax_ are assumed positive
 *  - @a n_ is a nonnegative Integer
 *  - @a b_ is normally positive, but may be negative.
 *  - Since @a b_ will generally be rather small, @a tmax_ should
 *    have a reasonable value (2 or 5) to make the routine
 *    more efficient.
 * @param[in] tmin_ Minimal allowed \f$t\f$
 * @param[in] tmax_ Maximal allowed \f$t\f$
 * @return Mandelstam variable \f$t\f$
 * @author Benno List
 * @date 22 Jan 1993
 * @date 15 Apr 1994
 * @date 17 Apr 2014
 */
double GenTDL(double tmin_, double tmax_, double b_, int n_);
/**
 * Generate the helicity of a photon
 * @param[in] longFr_ Fraction of longitydinally polarized photons
 * @return Helicity of the photon:
 *  - -1, +1: Transverse photon
 *  -  0: longitudinal photon
 * @author Benno List
 * @date 27 May 1993
 * @date 28 Apr 2014
 */
int Heli(double longFr_);
double ThetaToEta(double theta_);
double EtaToTheta(double eta_);
double EtaToY(double eta_, double m_, double pt_);

#endif
