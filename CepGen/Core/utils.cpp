#include "CepGen/Core/utils.h"
#include "CepGen/Core/Exception.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void Map( double expo, double xmin, double xmax, double& out, double& dout, const std::string& var_name_ )
{
  const double y = xmax/xmin;
  out = xmin*pow( y, expo );
  dout = out*log( y );
  DebuggingInsideLoop( "Map" )
    << "Mapping variable \"" << var_name_ << "\"\n\t"
    << "min = " << xmin << "\n\t"
    << "max = " << xmax << "\n\t"
    << "max/min = " << y << "\n\t"
    << "exponent = " << expo << "\n\t"
    << "output = " << out << "\n\t"
    << "d(output) = " << dout;
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

