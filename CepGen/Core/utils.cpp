#include "CepGen/Core/utils.h"
#include "CepGen/Core/Exception.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

