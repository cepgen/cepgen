#include "BlockDurandHa.h"

namespace CepGen
{
  namespace SF
  {
    BlockDurandHaParameters
    BlockDurandHaParameters::standard()
    {
      BlockDurandHaParameters p;
      p.a = { 8.205e-4, -5.148e-2, -4.725e-3 };
      p.b = { 2.214e-3,  1.244e-2,  5.958e-4 };
      p.c = { 0.255e0, 1.475e-1 };
      p.xn = 11.49;
      p.xlambda = 2.430;
      p.mu2 = 2.82;
      p.xm2 = 0.753;
      return p;
    }

    StructureFunctions BlockDurandHa( double q2, double xbj, const BlockDurandHaParameters& params )
    {
      StructureFunctions bdh;
      if ( q2 <= 0. ) return bdh;

      const double tau = q2 / ( q2 + params.mu2 );
      const double xl = log( 1. + q2 / params.mu2 );
      const double xlx = log( tau/xbj );

      const double AA = params.a[0] + params.a[1]*xl + params.a[2]*xl*xl;
      const double BB = params.b[0] + params.b[1]*xl + params.b[2]*xl*xl;
      const double CC = params.c[0] + params.c[1]*xl;
      const double DD = q2*( q2+params.xlambda*params.xm2 ) / pow( q2+params.xm2, 2 );

      bdh.F2 = DD*pow( 1.-xbj, params.xn ) * ( CC + AA*xlx + BB*xlx*xlx );

      return bdh;
    }
  }
}
