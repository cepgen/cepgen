#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Core/Exception.h"

#include <math.h>
#include <assert.h>
#include <vector>

namespace cepgen
{
  namespace strfun
  {
    /// \f$F_2\f$ parameterisation from Block, Durand, and Ha \cite Block:2014kza
    class BlockDurandHa : public Parameterisation
    {
      public:
        explicit BlockDurandHa( const ParametersList& params = ParametersList() );
        BlockDurandHa& operator()( double xbj, double q2 ) override;

      private:
        std::vector<double> a_, b_, c_;
        double n_;
        /// Effective mass spread parameter
        double lambda_;
        /// Asymptotic log-behaviour transition scale factor
        double mu2_;
        /// Squared effective mass (~VM mass)
        double m2_;
    };

    BlockDurandHa::BlockDurandHa( const ParametersList& params ) :
      Parameterisation( params ),
      a_( params.get<std::vector<double> >( "a", { 8.205e-4, -5.148e-2, -4.725e-3 } ) ),
      b_( params.get<std::vector<double> >( "b", { 2.217e-3,  1.244e-2,  5.958e-4 } ) ),
      c_( params.get<std::vector<double> >( "c", { 0.255e0, 1.475e-1 } ) ),
      n_     ( params.get<double>( "n", 11.49 ) ),
      lambda_( params.get<double>( "lambda", 2.430 ) ),
      mu2_   ( params.get<double>( "mu2", 2.82 ) ),
      m2_    ( params.get<double>( "m2", 0.753 ) )
    {
      assert( a_.size() == 3 );
      assert( b_.size() == 3 );
      assert( c_.size() == 2 );
    }

    BlockDurandHa&
    BlockDurandHa::operator()( double xbj, double q2 )
    {
      std::pair<double,double> nv = { xbj, q2 };
      if ( nv == old_vals_ )
        return *this;
      old_vals_ = nv;

      if ( q2 <= 0 ) {
        F2 = 0.;
        return *this;
      }

      const double tau = q2 / ( q2 + mu2_ );
      const double xl = log1p( q2 / mu2_ );
      const double xlx = log( tau/xbj );

      const double A = a_[0] + a_[1]*xl + a_[2]*xl*xl;
      const double B = b_[0] + b_[1]*xl + b_[2]*xl*xl;
      const double C = c_[0] + c_[1]*xl;
      const double D = q2*( q2+lambda_*m2_ ) / pow( q2+m2_, 2 );

      F2 = D*pow( 1.-xbj, n_ ) * ( C + A*xlx + B*xlx*xlx );

      return *this;
    }
  }
}

REGISTER_STRFUN( BlockDurandHa, strfun::BlockDurandHa )
