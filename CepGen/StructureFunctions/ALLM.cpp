#include "CepGen/StructureFunctions/ALLM.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include <cmath>
#include <cassert>

namespace cepgen
{
  namespace strfun
  {
    ALLM::ALLM( const ParametersList& params ) :
      Parameterisation( params ),
      params_( params.get<ParametersList>( "parameterisation" ) )
    {
      const auto& model = params.get<std::string>( "model" );
      if ( model == "GD07p" )
        params_ = Parameters::gd07p();
      else if ( model == "GD11p" )
        params_ = Parameters::gd11p();
      else if ( model == "ALLM91" )
        params_ = Parameters::allm91();
      else if ( model == "ALLM97" )
        params_ = Parameters::allm97();
      else if ( model == "HHT_ALLM" )
        params_ = Parameters::hht_allm();
      else if ( model == "HHT_ALLM_FT" )
        params_ = Parameters::hht_allm_ft();
      CG_DEBUG( "ALLM" ) << "ALLM structure functions builder initialised.\n"
        << "Parameterisation (" << params_.type << "):\n"
        << " *) Pomeron trajectory:\n"
        << "   a = {" << params_.pomeron.a.at( 0 ) << ", " << params_.pomeron.a.at( 1 ) << ", " << params_.pomeron.a.at( 2 ) << "}\n"
        << "   b = {" << params_.pomeron.b.at( 0 ) << ", " << params_.pomeron.b.at( 1 ) << ", " << params_.pomeron.b.at( 2 ) << "}\n"
        << "   c = {" << params_.pomeron.c.at( 0 ) << ", " << params_.pomeron.c.at( 1 ) << ", " << params_.pomeron.c.at( 2 ) << "}\n"
        << " *) Reggeon trajectory:\n"
        << "   a = {" << params_.reggeon.a.at( 0 ) << ", " << params_.reggeon.a.at( 1 ) << ", " << params_.reggeon.a.at( 2 ) << "}\n"
        << "   b = {" << params_.reggeon.b.at( 0 ) << ", " << params_.reggeon.b.at( 1 ) << ", " << params_.reggeon.b.at( 2 ) << "}\n"
        << "   c = {" << params_.reggeon.c.at( 0 ) << ", " << params_.reggeon.c.at( 1 ) << ", " << params_.reggeon.c.at( 2 ) << "}\n"
        << " masses: m₀²=" << params_.m02 << ", mp²=" << params_.mp2 << ", mr²=" << params_.mr2 << " GeV²\n"
        << " q₀²=" << params_.q02 << ", Λ²=" << params_.lambda2 << " GeV².";
    }

    ALLM&
    ALLM::operator()( double xbj, double q2 )
    {
      const double W2_eff = q2*( 1.-xbj )/xbj;
      const double xp = ( q2+params_.mp2 )/( q2+W2_eff+params_.mp2 ), xr = ( q2+params_.mr2 )/( q2+W2_eff+params_.mr2 );

      const double xlog1 = log( ( q2+params_.q02 )/ params_.lambda2 ), xlog2 = log( params_.q02/params_.lambda2 );
      const double t = log( xlog1/xlog2 );

      const double apom = params_.pomeron.a.at( 0 )+( params_.pomeron.a.at( 0 )-params_.pomeron.a.at( 1 ) )*( 1./( 1.+pow( t, params_.pomeron.a.at( 2 ) ) ) - 1. );
      const double bpom = params_.pomeron.b.at( 0 )+  params_.pomeron.b.at( 1 )*pow( t, params_.pomeron.b.at( 2 ) );
      const double cpom = params_.pomeron.c.at( 0 )+( params_.pomeron.c.at( 0 )-params_.pomeron.c.at( 1 ) )*( 1./( 1.+pow( t, params_.pomeron.c.at( 2 ) ) ) - 1. );

      const double areg = params_.reggeon.a.at( 0 )+params_.reggeon.a.at( 1 )*pow( t, params_.reggeon.a.at( 2 ) );
      const double breg = params_.reggeon.b.at( 0 )+params_.reggeon.b.at( 1 )*pow( t, params_.reggeon.b.at( 2 ) );
      const double creg = params_.reggeon.c.at( 0 )+params_.reggeon.c.at( 1 )*pow( t, params_.reggeon.c.at( 2 ) );

      const double F2_Pom = cpom*pow( xp, apom )*pow( 1.-xbj, bpom ), F2_Reg = creg*pow( xr, areg )*pow( 1.-xbj, breg );

      F2 = q2/( q2+params_.m02 ) * ( F2_Pom+F2_Reg );

      return *this;
    }

    //---------------------------------------------------------------------------------------------
    // parameterisation object
    //---------------------------------------------------------------------------------------------

    ALLM::Parameters::Parameters( const ParametersList& params ) :
      pomeron( params.get<ParametersList>( "pomeronTrajectory" ) ),
      reggeon( params.get<ParametersList>( "reggeonTrajectory" ) ),
      m02    ( params.get<double>( "m02" ) ),
      mp2    ( params.get<double>( "mp2" ) ),
      mr2    ( params.get<double>( "mr2" ) ),
      q02    ( params.get<double>( "q02" ) ),
      lambda2( params.get<double>( "lambda2" ) ),
      type   ( (Type)params.get<int>( "type", (int)Type::Invalid ) )
    {}

    ALLM::Parameters
    ALLM::Parameters::allm91()
    {
      static Parameters p( ParametersList()
        .set<ParametersList>( "pomeronTrajectory", ParametersList()
          .set<std::vector<double> >( "a", { -0.04503, -0.36407,  8.17091 } )
          .set<std::vector<double> >( "b", {  0.49222,  0.52116,  3.5515  } )
          .set<std::vector<double> >( "c", {  0.26550,  0.04856,  1.04682 } ) )
        .set<ParametersList>( "reggeonTrajectory", ParametersList()
          .set<std::vector<double> >( "a", {  0.60408,  0.17353,  1.61812 } )
          .set<std::vector<double> >( "b", {  1.26066,  1.83624,  0.81141 } )
          .set<std::vector<double> >( "c", {  0.67639,  0.49027,  2.66275 } ) )
        .set<double>( "m02", 0.30508 )
        .set<double>( "mp2", 10.676 )
        .set<double>( "mr2", 0.20623 )
        .set<double>( "q02", 0.27799 )
        .set<double>( "lambda2", 0.06527 )
        .set<int>( "type", (int)Type::ALLM91 ) );
      return p;
    }

    ALLM::Parameters
    ALLM::Parameters::allm97()
    {
      static Parameters p( ParametersList()
        .set<ParametersList>( "pomeronTrajectory", ParametersList()
          .set<std::vector<double> >( "a", { -0.0808,  -0.44812,  1.1709 } )
          .set<std::vector<double> >( "b", {  0.36292,  1.8917,   1.8439 } )
          .set<std::vector<double> >( "c", {  0.28067,  0.22291,  2.1979 } ) )
        .set<ParametersList>( "reggeonTrajectory", ParametersList()
          .set<std::vector<double> >( "a", {  0.58400,  0.37888,  2.6063  } )
          .set<std::vector<double> >( "b", {  0.01147,  3.7582,   0.49338 } )
          .set<std::vector<double> >( "c", {  0.80107,  0.97307,  3.4924  } ) )
        .set<double>( "m02", 0.31985 )
        .set<double>( "mp2", 49.457 )
        .set<double>( "mr2", 0.15052 )
        .set<double>( "q02", 0.52544 )
        .set<double>( "lambda2", 0.06526 )
        .set<int>( "type", (int)Type::ALLM97 ) );
      return p;
    }

    ALLM::Parameters
    ALLM::Parameters::hht_allm()
    {
      static Parameters p( ParametersList()
        .set<ParametersList>( "pomeronTrajectory", ParametersList()
          .set<std::vector<double> >( "a", { -0.835, -0.446,  10.6  } )
          .set<std::vector<double> >( "b", {-45.8,   55.7,    -0.031 } )
          .set<std::vector<double> >( "c", {  0.412,  0.164,  17.7  } ) )
        .set<ParametersList>( "reggeonTrajectory", ParametersList()
          .set<std::vector<double> >( "a", {  0.706,  0.185, -16.4  } )
          .set<std::vector<double> >( "b", { -1.29,   4.51,    1.16  } )
          .set<std::vector<double> >( "c", { -1.04,   2.97,    0.163 } ) )
        .set<double>( "m02", 0.446 )
        .set<double>( "mp2", 74.2 )
        .set<double>( "mr2", 29.3 )
        .set<double>( "q02", 4.74e-5 )
        .set<double>( "lambda2", 2.2e-8 ) );
      return p;
    }

    ALLM::Parameters
    ALLM::Parameters::hht_allm_ft()
    {
      static Parameters p( ParametersList()
        .set<ParametersList>( "pomeronTrajectory", ParametersList()
          .set<std::vector<double> >( "a", { -0.075, -0.470, 9.2   } )
          .set<std::vector<double> >( "b", { -0.477,  54.0,  0.073 } )
          .set<std::vector<double> >( "c", {  0.356,  0.171, 18.6  } ) )
        .set<ParametersList>( "reggeonTrajectory", ParametersList()
          .set<std::vector<double> >( "a", {  0.882, 0.082, -8.5   } )
          .set<std::vector<double> >( "b", {  0.339, 3.38,   1.07  } )
          .set<std::vector<double> >( "c", { -0.636, 3.37,  -0.660 } ) )
        .set<double>( "m02", 0.388 )
        .set<double>( "mp2", 50.8 )
        .set<double>( "mr2", 0.838 )
        .set<double>( "q02", 1.87e-5 )
        .set<double>( "lambda2", 4.4e-9 ) );
      return p;
    }

    ALLM::Parameters
    ALLM::Parameters::gd07p()
    {
      static Parameters p( ParametersList()
        .set<ParametersList>( "pomeronTrajectory", ParametersList()
          .set<std::vector<double> >( "a", { -0.105, -0.495, 1.29  } )
          .set<std::vector<double> >( "b", { -1.42,   4.51,  0.551 } )
          .set<std::vector<double> >( "c", {  0.339,  0.127, 1.16  } ) )
        .set<ParametersList>( "reggeonTrajectory", ParametersList()
          .set<std::vector<double> >( "a", { 0.374, 0.998, 0.775 } )
          .set<std::vector<double> >( "b", { 2.71,  1.83,  1.26  } )
          .set<std::vector<double> >( "c", { 0.838, 2.36,  1.77  } ) )
        .set<double>( "m02", 0.454 )
        .set<double>( "mp2", 30.7 )
        .set<double>( "mr2", 0.117 )
        .set<double>( "q02", 1.15 )
        .set<double>( "lambda2", 0.06527 )
        .set<int>( "type", (int)Type::GD07p ) );
      return p;
    }

    ALLM::Parameters
    ALLM::Parameters::gd11p()
    {
      static Parameters p( ParametersList()
        .set<ParametersList>( "pomeronTrajectory", ParametersList()
          .set<std::vector<double> >( "a", { -0.11895, -0.4783, 1.353 } )
          .set<std::vector<double> >( "b", {  1.0833,   2.656,  1.771 } )
          .set<std::vector<double> >( "c", {  0.3638,   0.1211, 1.166 } ) )
        .set<ParametersList>( "reggeonTrajectory", ParametersList()
          .set<std::vector<double> >( "a", {  0.3425,   1.0603, 0.5164 } )
          .set<std::vector<double> >( "b", {-10.408,   14.857,  0.07739 } )
          .set<std::vector<double> >( "c", {  1.3633,   2.256,  2.209   } ) )
        .set<double>( "m02", 0.5063 )
        .set<double>( "mp2", 34.75 )
        .set<double>( "mr2", 0.03190 )
        .set<double>( "q02", 1.374 )
        .set<double>( "lambda2", 0.06527 )
        .set<int>( "type", (int)Type::GD11p ) );
      return p;
    }

    ALLM::Parameters::Trajectory::Trajectory( const ParametersList& params ) :
      a( params.get<std::vector<double> >( "a", { 0., 0., 0. } ) ),
      b( params.get<std::vector<double> >( "b", { 0., 0., 0. } ) ),
      c( params.get<std::vector<double> >( "c", { 0., 0., 0. } ) )
    {
      assert( a.size() == 3 );
      assert( b.size() == 3 );
      assert( c.size() == 3 );
    }
  }
}

