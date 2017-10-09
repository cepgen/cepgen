#include "mstw_grid_handler.h"

namespace MSTW
{
  GridHandler::GridHandler( const char* filename ) :
    xacc_( 0 ), yacc_( 0 ), values_( 0 )
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

    values_ = new double[q2_vals.size()*xbj_vals.size()];
    for ( unsigned short i = 0; i < 2; ++i ) {
      splines_[i] = gsl_spline2d_alloc( T, q2_vals.size(), xbj_vals.size() );
    }
    xacc_ = gsl_interp_accel_alloc();
    yacc_ = gsl_interp_accel_alloc();

    std::cout << q2_vals.size() << "\t" << xbj_vals.size() << std::endl;
    while ( file >> q2 >> xbj >> f2 >> fl ) {
      gsl_spline2d_set( splines_[0], values_, q2, xbj, f2 );
      gsl_spline2d_set( splines_[1], values_, q2, xbj, fl );
    }
    std::vector<double> q2_vec( q2_vals.begin(), q2_vals.end() ), xbj_vec( xbj_vals.begin(), xbj_vals.end() );
    double* xa = &q2_vec[0], *ya = &xbj_vec[0];
    for ( unsigned short i = 0; i < 2; ++i ) {
      gsl_spline2d_init( splines_[i], xa, ya, values_, q2_vals.size(), xbj_vals.size() );
    }
  }

  GridHandler::~GridHandler()
  {
    for ( unsigned short i = 0; i < 2; ++i ) {
      gsl_spline2d_free( splines_[i] );
    }
    if ( xacc_ ) gsl_interp_accel_free( xacc_ );
    if ( yacc_ ) gsl_interp_accel_free( yacc_ );
    if ( values_ ) delete[] values_;
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
