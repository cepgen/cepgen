#include "BlockDurandHa.h"
#include <cmath>

namespace CepGen
{
  namespace SF
  {
    BlockDurandHa::Parameterisation
    BlockDurandHa::Parameterisation::standard()
    {
      Parameterisation p;
      p.a = { { 8.205e-4, -5.148e-2, -4.725e-3 } };
      p.b = { { 2.217e-3,  1.244e-2,  5.958e-4 } };
      p.c = { { 0.255e0, 1.475e-1 } };
      p.n = 11.49;
      p.lambda = 2.430;
      p.mu2 = 2.82;
      p.m2 = 0.753;
      return p;
    }

    BlockDurandHa::BlockDurandHa( const BlockDurandHa::Parameterisation& param ) :
      params_( param )
    {}

    BlockDurandHa&
    BlockDurandHa::operator()( double q2, double xbj )
    {
      std::pair<double,double> nv = { q2, xbj };
      if ( nv == old_vals_ )
        return *this;
      old_vals_ = nv;

      if ( q2 <= 0 ) {
        F2 = 0.;
        return *this;
      }

      const double tau = q2 / ( q2 + params_.mu2 );
      const double xl = log1p( q2 / params_.mu2 );
      const double xlx = log( tau/xbj );

      const double A = params_.a[0] + params_.a[1]*xl + params_.a[2]*xl*xl;
      const double B = params_.b[0] + params_.b[1]*xl + params_.b[2]*xl*xl;
      const double C = params_.c[0] + params_.c[1]*xl;
      const double D = q2*( q2+params_.lambda*params_.m2 ) / pow( q2+params_.m2, 2 );

      F2 = D*pow( 1.-xbj, params_.n ) * ( C + A*xlx + B*xlx*xlx );

      return *this;
    }
  }
}
