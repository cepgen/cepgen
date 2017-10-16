#include "CepGen/Core/utils.h"
#include "CepGen/Core/Exception.h"

#include <stdlib.h>
#include <stdarg.h>  // For va_start, etc.
#include <stdio.h>
#include <math.h>

void Map( double expo, double xmin, double xmax, double& out, double& dout, const std::string& var_name_ )
{
  const double y = xmax/xmin;
  out = xmin*pow( y, expo );
  dout = out*log( y );
  DebuggingInsideLoop( Form( "Mapping variable \"%s\"\n\t"
                             "min = %f\n\tmax = %f\n\tmax/min = %f\n\t"
                             "exponent = %f\n\t"
                             "output = %f\n\td(output) = %f",
                             var_name_.c_str(), xmin, xmax, y, expo, out, dout ) );
}

void Mapla( double y, double z, int u, double xm, double xp, double& x, double& d )
{
  double xmb, xpb, c, yy, zz, alp, alm, am, ap, ax;

  xmb = xm-y-z;
  xpb = xp-y-z;
  c = -4.*y*z;
  alp = sqrt( xpb*xpb + c );
  alm = sqrt( xmb*xmb + c );
  am = xmb+alm;
  ap = xpb+alp;
  yy = ap/am;
  zz = pow( yy, u );

  x = y + z + ( am*zz - c / ( am*zz ) ) / 2.;
  ax = sqrt( pow( x-y-z, 2 ) + c );
  d = ax*log( yy );
}

double BreitWigner( double er, double gamma, double emin, double emax, double x )
{
  if ( x==-1. ) x = drand();
  if ( gamma<1.e-3*er ) { return er; }

  const double a = atan( 2.*( emax-er ) / gamma ),
               b = atan( 2.*( emin-er ) / gamma ),
               e = er + gamma*tan( x*( a - b ) + b ) / 2.;

  if ( e>emax ) { return emax; }
  return e;
}

std::string
Form( const std::string fmt, ... )
{
  int size = ( (int)fmt.size() ) * 2 + 50;   // Use a rubric appropriate for your code
  std::string str;
  va_list ap;
  while ( true ) {     // Maximum two passes on a POSIX system...
    str.resize( size );
    va_start( ap, fmt );
    int n = vsnprintf( (char*)str.data(), size, fmt.c_str(), ap );
    va_end( ap );
    if ( n>-1 and n<size ) {  // Everything worked
      str.resize( n );
      return str;
    }
    if ( n>-1 )  // Needed size returned
      size = n + 1;   // For null char
    else size *= 2;      // Guess at a larger size (OS specific)
  }
  return str;
}

