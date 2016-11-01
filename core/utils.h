#ifndef utils_h
#define utils_h

#include <stdlib.h>
#include <stdarg.h>  // For va_start, etc.
#include <stdio.h>
#include <string.h>

#include "core/Exception.h"
#include "physics/Constants.h"

static std::stringstream error;

/// Provide a random number generated along a uniform distribution between 0 and 1
//inline double drand() { srand (time(NULL)); return static_cast<double>(rand())/RAND_MAX; }
#define drand() static_cast<double>( rand()/RAND_MAX )

/// Format a string using a printf style format descriptor.
std::string Form(const std::string fmt, ...);

inline const char* yesno( const bool& test ) { return ( test ) ? "\033[32;1myes\033[0m" : "\033[31;1mno\033[0m"; }
//inline const char* boldify( const char* str ) { const std::string out = std::string( "\033[33;1m" ) + std::string( str ) + std::string( "\033[0m" ); return out.c_str(); }
inline std::string boldify( const std::string& str ) { return Form( "\033[1m%s\033[0m", str.c_str() ); }
inline std::string boldify( const char* str ) { return boldify( std::string( str ) ); }
inline std::string boldify( const double& dbl ) { return boldify( Form("%.2f", dbl ) ); }
inline std::string boldify( const int& i ) { return boldify( Form("% d", i ) ); }
inline std::string boldify( const unsigned int& ui ) { return boldify( Form("%d", ui ) ); }
namespace Colour{
  enum TextColour { Gray=30, Red=31, Green=32, Yellow=33, Blue=34, Purple=35 };
}
inline std::string colourise( const std::string& str, const Colour::TextColour& col ) { return Form( "\033[%d%s\033[0m", col, str.c_str() ); }

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
double BreitWigner( double er, double gamma, double emin, double emax, double e=-1. );
/// Convert a polar angle to a pseudo-rapidity
inline double ThetaToEta( double theta_ ) { return -log( tan( theta_/180.*Constants::Pi/2. ) ); }
/// Convert a pseudo-rapidity to a polar angle
inline double EtaToTheta( double eta_ ) { return 2.*atan( exp( -eta_ ) )*180. / Constants::Pi; }
/// Convert a pseudo-rapidity to a rapidity
double EtaToY( double eta_, double m_, double pt_ );

#endif
