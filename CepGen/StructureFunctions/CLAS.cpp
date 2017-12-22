#include "CLAS.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  namespace SF
  {
    CLAS::Parameterisation
    CLAS::Parameterisation::standard()
    {
      Parameterisation params;
      params.mode = Parameterisation::proton;
      params.mp = ParticleProperties::mass( Proton );
      params.mpi0 = ParticleProperties::mass( PiZero );
      // SLAC fit parameters
      params.cn = { { 0.0640, 0.2250, 4.1060, -7.0790, 3.0550, 1.6421, 0.37636 } };
      params.cp = { { 0.25615, 2.1785, 0.89784, -6.7162, 3.7557, 1.6421, 0.37636 } };
      params.cd = { { 0.47709, 2.1602, 3.6274, -10.470, 4.9272, 1.5121, 0.35115 } };
      // CLAS parameterisation
      params.x = { { -0.599937, 4.76158, 0.411676 } };
      params.b = { { 0.755311, 3.35065, 3.51024, 1.74470 } };
      params.alpha = -0.174985;
      params.beta = 0.00967019;
      params.dmu = -0.0352567;
      params.dmup = 3.51852;
      params.ar = { { 1.04000, 0.481327, 0.655872, 0.747338 } };
      params.dmr = { { 1.22991, 1.51015, 1.71762, 1.95381 } };
      params.dgr = { { 0.106254, 0.0816620, 0.125520, 0.198915 } };
      params.lr = { { 1, 2, 3, 2 } };

      return params;
    }

    CLAS
    CLAS::operator()( double q2, double xbj ) const
    {
      const double mp2 = params_.mp*params_.mp;
      const double w2 = mp2 + q2*( 1.-xbj )/xbj;
      const double w_min = params_.mp+params_.mpi0;

      CLAS cl;
      if ( sqrt( w2 ) < w_min ) return cl;

      cl.F2 = f2slac( q2, xbj );

      double bkg = 0., resn = 0.;
      switch ( params_.mode ) {
        case Parameterisation::proton:
        case Parameterisation::neutron: {
          resbkg( q2, sqrt( w2 ), bkg, resn );
        } break;
        case Parameterisation::deuteron: {
        } break;
      }
      cl.F2 *= ( bkg+resn );
      return cl;
    }

    double
    CLAS::f2slac( double q2, double xbj ) const
    {
      if ( xbj >= 1. ) return 0.;
      double f2 = 0., xsxb = 0., xs = 0.;
      switch ( params_.mode ) {
        case Parameterisation::neutron: {
          xsxb = ( q2+params_.cn[6] )/( q2+params_.cn[5]*xbj );
          xs = xbj*xsxb;
          for ( unsigned short i = 0; i < 5; ++i ) {
            f2 += params_.cn[i]*pow( 1.-xs, i );
          }
        } break;
        case Parameterisation::proton: {
          xsxb = ( q2+params_.cp[6] )/( q2+params_.cp[5]*xbj );
          xs = xbj*xsxb;
          for ( unsigned short i = 0; i < 5; ++i ) {
            f2 += params_.cp[i]*pow( 1.-xs, i );
          }
        } break;
        case Parameterisation::deuteron: {
          xsxb = ( q2+params_.cd[6] )/( q2+params_.cd[5]*xbj );
          xs = xbj*xsxb;
          for ( unsigned short i = 0; i < 5; ++i ) {
            f2 += params_.cd[i]*pow( 1.-xs, i );
          }
          if ( xbj > 0. ) f2 /= ( 1.-exp( -7.70*( 1./xbj-1.+params_.mp*params_.mp/q2 ) ) );
        } break;
      }
      return f2 *= pow( 1.-xs, 3 ) / xsxb;
    }

    void
    CLAS::resbkg( double q2, double w, double& f2bkg, double& f2resn ) const
    {
      const double mp = params_.mp, mp2 = mp*mp;
      const double mpi0 = params_.mpi0, mpi02 = mpi0*mpi0;
      const double coef = 6.08974;
      // EVALUATE THE BACKGROUND AND RESONANCE TERMS OF THE
      // MODULATING FUNCTION (SLAC PARAMETRIZATION) FOR THE NUCLEON
      double wth = mp+mpi0;
      f2bkg = f2resn = 0.;
      if ( w < wth ) return;
      f2bkg = 1.;
      if ( w > 4. ) return;
      const double w2 = w*w;
      double qs = pow( w2+mp2-mpi02, 2 )-4.*mp2*w2;
      if ( qs <= 0. ) return;
      qs = 0.5 * sqrt( qs )/w;
      double omega = 0.5*( w2+q2-mp2 )/mp;
      const double xn = 0.5*q2/( mp*omega );
      const double bkg2 = ( w > params_.b[3] ) ? exp( -params_.b[2]*( w2-params_.b[3]*params_.b[3] ) ) : 1.;
      double etab = 1., etad = 1.;
      if ( q2 <= 2. && w <= 2.5 ) {
        etab = 1.-2.5*q2*exp( -12.5*q2*q2-50.*( w-1.325 )*( w-1.325 ) );
        etad = 1.+2.5*q2*exp( -12.5*q2*q2 );
      }
      f2bkg = params_.b[0]*( 1.-exp( -params_.b[1]*( w-wth ) ) )+( 1.-params_.b[0] )*( 1.-bkg2 );
      f2bkg *= etab*( 1.+( 1.-f2bkg )*( params_.x[0]+params_.x[1]*( xn-params_.x[2] )*( xn-params_.x[2] ) ) );
      double resn = 0.;
      for ( unsigned short i = 0; i < 4; ++i ) {
        double ai = params_.ar[i];
        if ( i == 0 ) ai = etad*( ai+q2*std::min( 0., params_.alpha+params_.beta*q2 ) );
        double dmi = params_.dmr[i];
        if ( i == 2 ) dmi *= ( 1.+params_.dmu/( 1.+params_.dmup*q2 ) );
        double qs0 = pow( dmi*dmi+mp2-mpi02, 2 )-4.*mp2*dmi*dmi;
        if ( qs0 <= 0. ) break;
        qs0 = 0.5*sqrt( qs0 )/dmi;
        int ji = 2*params_.lr[i];
        const double dg = 0.5*params_.dgr[i]*pow( qs/qs0, ji+1 )*( 1.+pow( coef*qs0, ji ) )/( 1.+pow( coef*qs, ji ) );
        resn += ai*dg/( ( w-dmi )*( w-dmi )+dg*dg );
      }
      resn *= 0.5*( 1.-params_.b[0] )*bkg2/( mp*M_PI );
    }
  }
}

