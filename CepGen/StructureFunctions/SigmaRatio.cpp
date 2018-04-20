#include "SigmaRatio.h"
#include "CepGen/Physics/ParticleProperties.h"
#include <math.h>
#include <iostream>

namespace CepGen
{
  namespace SF
  {
    const double SigmaRatio::mp_ = ParticleProperties::mass( PDG::Proton );
    const double SigmaRatio::mp2_ = SigmaRatio::mp_*SigmaRatio::mp_;

    double
    SigmaRatio::theta( double q2, double xbj ) const
    {
      return 1.+12.*( q2/( q2+1. ) )*( 0.125*0.125/( 0.125*0.125+xbj*xbj ) );
    }

    E143Ratio::Parameterisation
    E143Ratio::Parameterisation::standard()
    {
      Parameterisation out;
      out.q2_b = 0.34;
      out.lambda2 = 0.2*0.2;
      out.a = { { 0.0485, 0.5470,  2.0621, -0.3804,   0.5090, -0.0285 } };
      out.b = { { 0.0481, 0.6114, -0.3509, -0.4611,   0.7172, -0.0317 } };
      out.c = { { 0.0577, 0.4644,  1.8288, 12.3708, -43.1043, 41.7415 } };
      return out;
    }

    E143Ratio::E143Ratio( const Parameterisation& param ) :
      params_( param )
    {}

    double
    E143Ratio::operator()( double q2, double xbj, double& err ) const
    {
      const double u = q2/params_.q2_b;
      const double inv_xl = 1./log( q2/params_.lambda2 );
      const double pa = ( 1.+params_.a[3]*xbj+params_.a[4]*xbj*xbj )*pow( xbj, params_.a[5] );
      const double pb = ( 1.+params_.b[3]*xbj+params_.b[4]*xbj*xbj )*pow( xbj, params_.b[5] );
      const double theta = SigmaRatio::theta( q2, xbj );
      const double q2_thr = params_.c[3]*xbj + params_.c[4]*xbj*xbj+params_.c[5]*xbj*xbj*xbj;
      // here come the three fits
      const double ra = params_.a[0]*inv_xl*theta + params_.a[1]/pow( pow( q2, 4 )+pow( params_.a[2], 4 ), 0.25 )*pa,
                   rb = params_.b[0]*inv_xl*theta + ( params_.b[1]/q2+params_.b[2]/( q2*q2+0.3*0.3 ) )*pb,
                   rc = params_.c[0]*inv_xl*theta + params_.c[1]*pow( pow( q2-q2_thr, 2 )+pow( params_.c[2], 2 ), -0.5 );

      const double r = ( ra+rb+rc ) / 3.; // R is set to be the average of the three fits
      // numerical safety for low-Q²
      err = 0.0078-0.013*xbj+( 0.070-0.39*xbj+0.70*xbj*xbj )/( 1.7+q2 );
      if ( q2 > params_.q2_b ) return r;
      return r * 0.5 * ( 3.*u-u*u*u );
    }

    R1990Ratio::Parameterisation
    R1990Ratio::Parameterisation::standard()
    {
      Parameterisation out;
      out.lambda2 = 0.2*0.2;
      out.b = { { 0.0635, 0.5747, -0.3534 } };
      return out;
    }

    R1990Ratio::R1990Ratio( const Parameterisation& param ) :
      params_( param )
    {}

    double
    R1990Ratio::operator()( double q2, double xbj, double& err ) const
    {
      err = 0.;
      return ( params_.b[0]+SigmaRatio::theta( q2, xbj )/log( q2/params_.lambda2 )
             + params_.b[1]/q2
             + params_.b[2]/( q2*q2+0.09 ) );
    }

    double
    CLASRatio::operator()( double q2, double xbj, double& err ) const
    {
      // 2 kinematic regions:
      //  - resonances ( w < 2.5 )
      //  - DIS ( w > 2.5 )
      const double w2 = mp2_ + q2*( 1.-xbj )/xbj, w = sqrt( w2 );
      const double xth = q2/( q2+2.5*2.5-mp2_ ); // xth = x( W = 2.5 GeV )
      const double zeta = log( 25.*q2 );
      const double xitmp = ( w < 2.5 ) ? theta( q2, xth ) : theta( q2, xbj );
      const double tmp = 0.041*xitmp/zeta + 0.592/q2 - 0.331/( 0.09+q2*q2 );
      if ( w < 2.5 ) return tmp * pow( ( 1.-xbj )/( 1.-xth ), 3 );
      return tmp;
    }

    double
    SBRatio::operator()( double q2, double xbj, double& err ) const
    {
      err = 0.;
      return 0.014*q2*( exp( -0.07*q2 )+41.*exp( -0.8*q2 ) );
    }
  }
}

