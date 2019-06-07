#include "CepGen/StructureFunctions/SigmaRatio.h"

#include "CepGen/Physics/PDG.h"

#include <iostream>
#include <cmath>
#include <cassert>

namespace cepgen
{
  namespace sigrat
  {
    const double Parameterisation::mp_ = PDGInfo::get()( PDG::proton ).mass;
    const double Parameterisation::mp2_ = Parameterisation::mp_*Parameterisation::mp_;

    Parameterisation::Parameterisation( const ParametersList& params )
    {}

    double
    Parameterisation::theta( double xbj, double q2 ) const
    {
      return 1.+12.*( q2/( q2+1. ) )*( 0.125*0.125/( 0.125*0.125+xbj*xbj ) );
    }

    //---------------------------------------------------------------------------------------------

    E143::E143( const ParametersList& params ) :
      q2_b_   ( params.get<double>( "q2_b", 0.34 ) ),
      lambda2_( params.get<double>( "lambda2", 0.2*0.2 ) ),
      a_      ( params.get<std::vector<double> >( "a", { 0.0485, 0.5470,  2.0621, -0.3804,   0.5090, -0.0285 } ) ),
      b_      ( params.get<std::vector<double> >( "b", { 0.0481, 0.6114, -0.3509, -0.4611,   0.7172, -0.0317 } ) ),
      c_      ( params.get<std::vector<double> >( "c", { 0.0577, 0.4644,  1.8288, 12.3708, -43.1043, 41.7415 } ) )
    {
      assert( a_.size() == 6 );
      assert( b_.size() == 6 );
      assert( c_.size() == 6 );
    }

    double
    E143::operator()( double xbj, double q2, double& err ) const
    {
      const double u = q2/q2_b_;
      const double inv_xl = 1./log( q2/lambda2_ );
      const double pa = ( 1.+a_.at( 3 )*xbj+a_.at( 4 )*xbj*xbj )*pow( xbj, a_.at( 5 ) );
      const double pb = ( 1.+b_.at( 3 )*xbj+b_.at( 4 )*xbj*xbj )*pow( xbj, b_.at( 5 ) );
      const double q2_thr =  c_.at( 3 )*xbj+c_.at( 4 )*xbj*xbj+c_.at( 5 )*xbj*xbj*xbj;
      const double th = theta( xbj, q2 );
      // here come the three fits
      const double ra = a_.at( 0 )*inv_xl*th+  a_.at( 1 )/pow( pow( q2, 4 )+pow( a_.at( 2 ), 4 ), 0.25 )*pa,
                   rb = b_.at( 0 )*inv_xl*th+( b_.at( 1 )/q2+b_.at( 2 )/( q2*q2+0.3*0.3 ) )*pb,
                   rc = c_.at( 0 )*inv_xl*th+  c_.at( 1 )/std::hypot( q2-q2_thr, c_.at( 2 ) );

      const double r = ( ra+rb+rc ) / 3.; // R is set to be the average of the three fits
      // numerical safety for low-Q²
      err = 0.0078-0.013*xbj+( 0.070-0.39*xbj+0.70*xbj*xbj )/( 1.7+q2 );
      if ( q2 > q2_b_ )
        return r;
      return r * 0.5 * ( 3.*u-u*u*u );
    }

    //---------------------------------------------------------------------------------------------

    R1990::R1990( const ParametersList& params ) :
      lambda2_( params.get<double>( "lambda2", 0.04 ) ),
      b_      ( params.get<std::vector<double> >( "b", { 0.0635, 0.5747, -0.3534 } ) )
    {
      assert( b_.size() == 3 );
    }

    double
    R1990::operator()( double xbj, double q2, double& err ) const
    {
      err = 0.;
      return b_.at( 0 )+theta( xbj, q2 )/log( q2/lambda2_ )+b_.at( 1 )/q2+b_.at( 2 )/( q2*q2+0.09 );
    }

    //---------------------------------------------------------------------------------------------

    CLAS::CLAS( const ParametersList& params ) :
      p_  ( params.get<std::vector<double> >( "p", { 0.041, 0.592, 0.331 } ) ),
      wth_( params.get<double>( "wth", 2.5 ) ),
      q20_( params.get<double>( "q20", 0.3 ) )
    {
      assert( p_.size() == 3 );
    }

    double
    CLAS::operator()( double xbj, double q2, double& err ) const
    {
      //--- 2 kinematic regions: resonances ( w < wth ), and DIS ( w > wth )
      const double w2 = mp2_ + q2*( 1.-xbj )/xbj, w = sqrt( w2 );
      const double xth = q2/( q2+wth_*wth_-mp2_ ); // xth = x( W = wth )
      const double zeta = log( 25.*q2 );
      const double xitmp = ( w < wth_ ) ? theta( xth, q2 ) : theta( xbj, q2 );
      const double tmp = p_.at( 0 )*xitmp/zeta + p_.at( 1 )/q2 - p_.at( 3 )/( q20_*q20_+q2*q2 );
      if ( w >= wth_ )
        return tmp;
      return tmp * pow( ( 1.-xbj )/( 1.-xth ), 3 );
    }

    //---------------------------------------------------------------------------------------------

    SibirtsevBlunden::SibirtsevBlunden( const ParametersList& params ) :
      a_ ( params.get<double>( "a", 0.014 ) ),
      b1_( params.get<double>( "b1", -0.07 ) ),
      b2_( params.get<double>( "b2", -0.8 ) ),
      c_ ( params.get<double>( "c", 41. ) )
    {}

    double
    SibirtsevBlunden::operator()( double xbj, double q2, double& err ) const
    {
      err = 0.;
      //--- equation (10) of reference paper
      return a_*q2*( exp( b1_*q2 )+c_*exp( b2_*q2 ) );
    }
  }
}

