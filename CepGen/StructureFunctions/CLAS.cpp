#include "CepGen/StructureFunctions/CLAS.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Constants.h"

#include "CepGen/Event/Particle.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  namespace strfun
  {
    CLAS::Parameters
    CLAS::Parameters::standard_proton()
    {
      Parameters params;
      params.mode = Parameters::proton;
      params.mp = mp_;
      params.mpi0 = particleproperties::mass( PDG::piZero );
      // SLAC fit parameters
      params.c_slac = { { 0.25615, 2.1785, 0.89784, -6.7162, 3.7557, 1.6421, 0.37636 } };
      // CLAS parameterisation
      params.x = { { -0.599937, 4.76158, 0.411676 } };
      params.b = { { 0.755311, 3.35065, 3.51024, 1.74470 } };
      params.alpha = -0.174985;
      params.beta = 0.00967019;
      params.mu = -0.0352567;
      params.mup = 3.51852;

      Parameters::Resonance r0;
      r0.amplitude = 1.04;
      r0.mass = 1.22991;
      r0.width = 0.106254;
      r0.angular_momentum = 1;
      params.resonances.emplace_back( r0 );

      Parameters::Resonance r1;
      r1.amplitude = 0.481327;
      r1.mass = 1.51015;
      r1.width = 0.0816620;
      r1.angular_momentum = 2;
      params.resonances.emplace_back( r1 );

      Parameters::Resonance r2;
      r2.amplitude = 0.655872;
      r2.mass = 1.71762;
      r2.width = 0.125520;
      r2.angular_momentum = 3;
      params.resonances.emplace_back( r2 );

      Parameters::Resonance r3;
      r3.amplitude = 0.747338;
      r3.mass = 1.95381;
      r3.width = 0.198915;
      r3.angular_momentum = 2;
      params.resonances.emplace_back( r3 );

      return params;
    }

    CLAS::Parameters
    CLAS::Parameters::standard_neutron()
    {
      Parameters params = standard_proton();
      params.mode = Parameters::neutron;
      params.c_slac = { { 0.0640, 0.2250, 4.1060, -7.0790, 3.0550, 1.6421, 0.37636 } };
      return params;
    }

    CLAS::Parameters
    CLAS::Parameters::standard_deuteron()
    {
      Parameters params = standard_proton();
      params.mode = Parameters::deuteron;
      params.c_slac = { { 0.47709, 2.1602, 3.6274, -10.470, 4.9272, 1.5121, 0.35115 } };
      params.x = { { -0.21262, 6.9690, 0.40314 } };
      params.b = { { 0.76111, 4.1470, 3.7119, 1.4218 } };
      params.alpha = -0.24480;
      params.beta = 0.014503;

      params.resonances.clear();

      Parameters::Resonance r0;
      r0.amplitude = 0.74847;
      r0.mass = 1.2400;
      r0.width = 0.12115;
      r0.angular_momentum = 1;
      params.resonances.emplace_back( r0 );

      Parameters::Resonance r1;
      r1.amplitude = 0.011500;
      r1.mass = 1.4772;
      r1.width = 0.0069580;
      r1.angular_momentum = 2;
      params.resonances.emplace_back( r1 );

      Parameters::Resonance r2;
      r2.amplitude = 0.12662;
      r2.mass = 1.5233;
      r2.width = 0.084095;
      r2.angular_momentum = 3;
      params.resonances.emplace_back( r2 );

      Parameters::Resonance r3;
      r3.amplitude = 0.747338;
      r3.mass = 1.95381;
      r3.width = 0.198915;
      r3.angular_momentum = 2;
      params.resonances.emplace_back( r3 );

      return params;
    }

    CLAS::CLAS( const Parameters& params ) :
      Parameterisation( Type::CLAS ), params_( params )
    {}

    CLAS&
    CLAS::operator()( double xbj, double q2 )
    {
      std::pair<double,double> nv = { xbj, q2 };
      if ( nv == old_vals_ )
        return *this;
      old_vals_ = nv;

      const double mp2 = params_.mp*params_.mp;
      const double w2 = mp2 + q2*( 1.-xbj )/xbj;
      const double w_min = params_.mp+params_.mpi0;

      if ( sqrt( w2 ) < w_min ) {
        F2 = 0.;
        return *this;
      }

      F2 = f2slac( xbj, q2 );
      std::pair<double,double> rb = resbkg( q2, sqrt( w2 ) );

      F2 *= ( rb.first+rb.second );
      return *this;
    }

    double
    CLAS::f2slac( double xbj, double q2 ) const
    {
      if ( xbj >= 1. )
        return 0.;

      const double xsxb = ( q2+params_.c_slac[6] )/( q2+params_.c_slac[5]*xbj );
      const double xs = xbj*xsxb;

      double f2 = 0.;
      for ( unsigned short i = 0; i < 5; ++i )
        f2 += params_.c_slac[i]*pow( 1.-xs, i );

      if ( params_.mode == Parameters::deuteron && xbj > 0. )
        f2 /= ( 1.-exp( -7.70*( 1./xbj-1.+params_.mp*params_.mp/q2 ) ) );

      return f2 * pow( 1.-xs, 3 ) / xsxb;
    }

    std::pair<double,double>
    CLAS::resbkg( double q2, double w ) const
    {
      const double mp2 = params_.mp*params_.mp, mpi02 = params_.mpi0*params_.mpi0;
      const double coef = 6.08974;

      double wth = params_.mp+params_.mpi0;
      if ( w < wth )
        return std::make_pair( 0., 0. );
      if ( w > 4. )
        return std::make_pair( 1., 0. );

      const double w2 = w*w;

      double qs = pow( w2+mp2-mpi02, 2 )-4.*mp2*w2;
      if ( qs <= 0. ) return std::make_pair( 1., 0. );
      qs = 0.5 * sqrt( qs )/w;

      const double omega = 0.5*( w2+q2-mp2 )/params_.mp;
      const double xn = 0.5*q2/( params_.mp*omega );

      const double bkg2 = ( w > params_.b[3] )
        ? exp( -params_.b[2]*( w2-params_.b[3]*params_.b[3] ) )
        : 1.;

      double f2bkg = (    params_.b[0] )*( 1.-exp( -params_.b[1]*( w-wth ) ) )
                   + ( 1.-params_.b[0] )*( 1.-bkg2 );
      f2bkg *= ( 1.+( 1.-f2bkg )*( params_.x[0]+params_.x[1]*pow( xn-params_.x[2], 2 ) ) );

      double etab = 1., etad = 1.;
      if ( params_.mode != Parameters::deuteron && q2 <= 2. && w <= 2.5 ) {
        etab = 1.-2.5*q2*exp( -12.5*q2*q2-50.*( w-1.325 )*( w-1.325 ) );
        etad = 1.+2.5*q2*exp( -12.5*q2*q2 );
      }
      f2bkg *= etab;

      double f2resn = 0.;

      for ( unsigned short i = 0; i < params_.resonances.size(); ++i ) {
        const Parameters::Resonance& res = params_.resonances[i];
        const double ai = ( i == 0 )
          ? etad * ( res.amplitude + q2*std::min( 0., params_.alpha+params_.beta*q2 ) )
          : res.amplitude;
        const double dmi = ( i == 2 )
          ? res.mass * ( 1.+params_.mu/( 1.+params_.mup*q2 ) )
          : res.mass;
        double qs0 = pow( dmi*dmi+mp2-mpi02, 2 )-4.*mp2*dmi*dmi;
        if ( qs0 <= 0. )
          break;
        qs0 = 0.5*sqrt( qs0 )/dmi;
        int ji = 2*res.angular_momentum;
        const double dg = 0.5*res.width*pow( qs/qs0, ji+1 )*( 1.+pow( coef*qs0, ji ) )/( 1.+pow( coef*qs, ji ) );
        f2resn += ai*dg/( ( w-dmi )*( w-dmi )+dg*dg );
      }
      f2resn *= 0.5*( 1.-params_.b[0] )*bkg2/params_.mp*M_1_PI;

      return std::make_pair( f2bkg, f2resn );
    }
  }
}
