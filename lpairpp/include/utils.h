#ifndef _UTILS_H
#define _UTILS_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <fstream>

/** @brief Electromagnetic coupling constant \f$\alpha_{em}=\frac{e^2}{4\pi\epsilon_0\hbar c}\f$ */
#define alphaF 1./137.04
/** @brief \f$\frac{1}{(\hbar c)^2}~[\mathrm b^{-1}]\f$? */
#define muBarn 1./389.39
#define pi 3.1415926535897932384626434
#define sconst 3.89351824E8
#define sconstb 2.1868465E10

/**
 * Gets the mass in GeV/c**2 of a particle given its PDG identifier
 * @brief Gets the mass of a particle
 * @param pdgId_ PDG ID of the particle whose mass is requested
 * @return Mass of the particle in GeV/c**2
 */
double GetMassFromPDGId(int pdgId_);

/**
 * Computes the proton structure function (F.W Brasse et al., DESY 76/11 (1976),
 *   http://dx.doi.org/10.1016/0550-3213(76)90231-5)
 * @cite Brasse1976413
 */
bool PSF(double,double,double*,double*,double*);
/**
 * @brief Defines modified variables of integration to avoid peaks integrations (see @cite Vermaseren1983347 for details)
 * Returns a set of two modified variables of integration to maintain the stability of the integrant. These two new variables are :
 * - \f$y_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\f$ the new variable
 * - \f$\mathrm dy_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\log\frac{x_{min}}{x_{max}}\f$, the new variable's differential form
 * @brief Redefines the variables of integration in order to avoid the strong peaking of the integrant.
 * @param expo_ Exponant
 * @param xmin_ Minimal value of the variable
 * @param xmax_ Maximal value of the variable
 * @param out_ The new variable definition
 * @param dout_ The differential variant of the new variable definition
 * @note This method overrides the set of `mapxx` subroutines in ILPAIR, with a slight difference according to the sign of the
 *  \f$\mathrm dy_{out}\f$ parameter :
 *  - left unchanged :
 * > `mapw2`, `mapxq`, `mapwx`, `maps2`
 *  - opposite sign :
 * > `mapt1`, `mapt2`
 */
void Map(double expo_, double xmin_, double xmax_, double* out_, double* dout_);
void Mapla(double,double,int,double,double,double*,double*);
//void Symmetrise(double, double, double*, double*);
void Lorenb(double u_, double ps_[], double pi_[], double pf_[]);

/**
 * Generate random number with Breit-Wigner distribution
 * @return Random number between emin_ and emax_ with Breit-Wigner distribution:
 *  \f$\frac{1}{(E-E_r)^2+\Gamma^2/4}\f$
 * @param er_ Maximum of distribution
 * @param gamma_ Width of distribution
 * @param emin_ Minimal value of RanBW
 * @param emax_ Maximal value of RanBW
 */
double RanBW(double er_, double gamma_, double emin_, double emax_);
double GenerT(double tmin_, double tmax_, double b_, double anexp_);

#endif
