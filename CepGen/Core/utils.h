#ifndef CepGen_Core_utils_h
#define CepGen_Core_utils_h

#include <stdlib.h>
#include <stdarg.h>  // For va_start, etc.
#include <stdio.h>
#include <string.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/Constants.h"

static std::stringstream error;

/// Provide a random number generated along a uniform distribution between 0 and 1
//inline double drand() { srand (time(nullptr)); return static_cast<double>(rand())/RAND_MAX; }
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
 * @param[in] expo Exponant
 * @param[in] xmin Minimal value of the variable
 * @param[in] xmax Maximal value of the variable
 * @param[out] out The new variable definition
 * @param[out] dout The differential variant of the new variable definition
 * @param[in] var_name The variable name
 * @note This method overrides the set of `mapxx` subroutines in ILPAIR, with a slight difference according to the sign of the
 *  \f$\mathrm dy_{out}\f$ parameter :
 *  - left unchanged :
 * > `mapw2`, `mapxq`, `mapwx`, `maps2`
 *  - opposite sign :
 * > `mapt1`, `mapt2`
 */
void Map( double expo, double xmin, double xmax, double& out, double& dout, const std::string& var_name="" );
void Mapla( double y, double z, int u, double xm, double xp, double& x, double& d );

#endif
