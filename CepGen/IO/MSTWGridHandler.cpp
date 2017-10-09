#include "MSTWGridHandler.h"

namespace MSTW
{
  GridHandler::GridHandler( const char* filename ) :
    xacc_( 0 ), yacc_( 0 ), values_( { 0, 0 } )
  {
    const gsl_interp2d_type* T = gsl_interp2d_bilinear;
    std::ifstream file( filename );
    if ( !file.is_open() ) {
      FatalError( Form( "Impossible to load grid file \"%s\"!", filename ) );
    }

    double q2, xbj, f2, fl;

    // first loop to evaluate the limits
    std::set<double> q2_vals, xbj_vals;
    while ( file >> q2 >> xbj >> f2 >> fl ) {
      q2_vals.insert( q2 );
      xbj_vals.insert( xbj );
    }

    if ( q2_vals.size() < 2 || xbj_vals.size() < 2 ) {
      FatalError( "Invalid grid retrieved!" );
    }
    Information( Form( "Q2 in range [%.3e:%.3e]\n\txBj in range [%.3e:%.3e]", *q2_vals.begin(), *q2_vals.rbegin(), *xbj_vals.begin(), *xbj_vals.rbegin() ) );

    for ( unsigned short i = 0; i < 2; ++i ) {
      values_[i] = new double[q2_vals.size()*xbj_vals.size()];
      splines_[i] = gsl_spline2d_alloc( T, q2_vals.size(), xbj_vals.size() );
    }
    xacc_ = gsl_interp_accel_alloc();
    yacc_ = gsl_interp_accel_alloc();

    // second loop to populate the grid
    file.clear(); file.seekg( 0, std::ios::beg );
    while ( file >> q2 >> xbj >> f2 >> fl ) {
      unsigned short id_q2 = std::distance( q2_vals.begin(), q2_vals.lower_bound( q2 ) ),
                     id_xbj = std::distance( xbj_vals.begin(), xbj_vals.lower_bound( xbj ) );
      gsl_spline2d_set( splines_[0], values_[0], id_q2, id_xbj, f2 );
      gsl_spline2d_set( splines_[1], values_[1], id_q2, id_xbj, fl );
    }

    // initialise the splines object
    std::vector<double> q2_vec( q2_vals.begin(), q2_vals.end() ), xbj_vec( xbj_vals.begin(), xbj_vals.end() );
    double* xa = &q2_vec[0], *ya = &xbj_vec[0];
    for ( unsigned short i = 0; i < 2; ++i ) {
      gsl_spline2d_init( splines_[i], xa, ya, values_[i], q2_vals.size(), xbj_vals.size() );
    }
  }

  GridHandler::~GridHandler()
  {
    for ( unsigned short i = 0; i < 2; ++i ) {
      gsl_spline2d_free( splines_[i] );
      if ( values_[i] ) delete[] values_[i];
    }
    if ( xacc_ ) gsl_interp_accel_free( xacc_ );
    if ( yacc_ ) gsl_interp_accel_free( yacc_ );
  }

  CepGen::StructureFunctions
  GridHandler::eval( double q2, double xbj ) const
  {
    CepGen::StructureFunctions ev;
    ev.F2 = gsl_spline2d_eval( splines_[0], q2, xbj, xacc_, yacc_ );
    ev.FL = gsl_spline2d_eval( splines_[1], q2, xbj, xacc_, yacc_ );
    return ev;
  }
}
