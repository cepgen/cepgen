#include "SigmaRatio.h"
#include "CepGen/Physics/ParticleProperties.h"
#include <math.h>
#include <iostream>

namespace CepGen
{
  namespace SF
  {
    E143Ratio::Parameterisation
    E143Ratio::Parameterisation::standard()
    {
      Parameterisation out;
      out.q2_b = 0.34;
      out.a = { { 0.0485, 0.5470,  2.0621, -0.3804,   0.5090, -0.0285 } };
      out.b = { { 0.0481, 0.6114, -0.3509, -0.4611,   0.7172, -0.0317 } };
      out.c = { { 0.0577, 0.4644,  1.8288, 12.3708, -43.1043, 41.7415 } };
      return out;
    }

    double
    E143Ratio::operator()( double q2, double xbj ) const
    {
      const double u = q2/params_.q2_b;
      const double xl = log( 25.*q2 );
      const double pa = ( 1.+params_.a[3]*xbj+params_.a[4]*xbj*xbj )*pow( xbj, params_.a[5] );
      const double pb = ( 1.+params_.b[3]*xbj+params_.b[4]*xbj*xbj )*pow( xbj, params_.b[5] );
      const double tt = xi( q2, xbj );
      const double q2_thr = params_.c[3]*xbj + params_.c[4]*xbj*xbj+params_.c[5]*xbj*xbj*xbj;
      const double ra = params_.a[0]/xl*tt + params_.a[1]/pow( pow( q2, 4 )+pow( params_.a[2], 4 ), 0.25 )*pa,
                   rb = params_.b[0]/xl*tt + ( params_.b[1]/q2+params_.b[2]/( q2*q2+0.3*0.3 ) )*pb,
                   rc = params_.c[0]/xl*tt + params_.c[1]*pow( pow( q2-q2_thr, 2 )+pow( params_.c[2], 2 ), -0.5 );

      const double r = ( ra+rb+rc ) / 3.;
      // numerical safety for low-Q²
      if ( q2 > params_.q2_b ) return r;
      return r * 0.5 * ( 3.*u-u*u*u );
    }

    double
    CLASRatio::operator()( double q2, double xbj ) const
    {
      // 2 kinematic regions:
      //  - resonances ( w < 2.5 )
      //  - DIS ( w > 2.5 )
      const double mp = ParticleProperties::mass( Proton );
      const double w2 = mp*mp + q2*( 1.-xbj )/xbj, w = sqrt( w2 );
      const double xth = q2/( q2+2.5*2.5-mp*mp ); // xth = x( W = 2.5 GeV )
      const double zeta = log( 25.*q2 );
      const double xitmp = ( w < 2.5 ) ? xi( q2, xth ) : xi( q2, xbj );
      const double tmp = 0.041*xitmp/zeta + 0.592/q2 - 0.331/( 0.09+q2*q2 );
      if ( w < 2.5 ) return tmp * pow( ( 1.-xbj )/( 1.-xth ), 3 );
      return tmp;
    }

    double
    SigmaRatio::xi( double q2, double xbj ) const
    {
      return 1.+12.*( q2/( q2+1. ) )*( 0.125*0.125/( 0.125*0.125+xbj*xbj ) );
    }
  }
}

