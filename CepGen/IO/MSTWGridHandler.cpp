#include "MSTWGridHandler.h"

namespace MSTW
{
  GridHandler&
  GridHandler::get( const char* filename )
  {
    static GridHandler instance( filename );
    return instance;
  }

  GridHandler::GridHandler( const char* filename )
#ifdef GOOD_GSL
    : xacc_( 0 ), yacc_( 0 ), values_( { 0, 0 } )
#endif
  {
    gsl_set_error_handler_off();
#ifdef GOOD_GSL
    const gsl_interp2d_type* T = gsl_interp2d_bilinear;
    std::ifstream file( filename, std::ios::binary | std::ios::in );
    if ( !file.is_open() ) {
      FatalError( Form( "Impossible to load grid file \"%s\"!", filename ) );
    }

    sfval_t val;

    // first loop to evaluate the limits
    std::set<double> q2_vals, xbj_vals;
    while ( file.read( reinterpret_cast<char*>( &val ), sizeof( sfval_t ) ) ) {
      q2_vals.insert( val.q2 );
      xbj_vals.insert( val.xbj );
    }

    if ( q2_vals.size() < 2 || xbj_vals.size() < 2 ) {
      FatalError( "Invalid grid retrieved!" );
    }
    Information( Form( "MSTW structure functions grid evaluator built\n\t"
                       " Q² in range [%.3e:%.3e]\n\t"
                       "xBj in range [%.3e:%.3e]",
                       *q2_vals.begin(), *q2_vals.rbegin(), *xbj_vals.begin(), *xbj_vals.rbegin() ) );

    for ( unsigned short i = 0; i < 2; ++i ) {
      values_[i] = new double[q2_vals.size()*xbj_vals.size()];
      splines_[i] = gsl_spline2d_alloc( T, q2_vals.size(), xbj_vals.size() );
    }
    xacc_ = gsl_interp_accel_alloc();
    yacc_ = gsl_interp_accel_alloc();

    // second loop to populate the grid
    file.clear(); file.seekg( 0, std::ios::beg );
    while ( file.read( reinterpret_cast<char*>( &val ), sizeof( sfval_t ) ) ) {
      unsigned short id_q2 = std::distance( q2_vals.begin(), q2_vals.lower_bound( val.q2 ) ),
                     id_xbj = std::distance( xbj_vals.begin(), xbj_vals.lower_bound( val.xbj ) );
      gsl_spline2d_set( splines_[0], values_[0], id_q2, id_xbj, val.f2 );
      gsl_spline2d_set( splines_[1], values_[1], id_q2, id_xbj, val.fl );
    }

    // initialise the splines object
    std::vector<double> q2_vec( q2_vals.begin(), q2_vals.end() ), xbj_vec( xbj_vals.begin(), xbj_vals.end() );
    double* xa = &q2_vec[0], *ya = &xbj_vec[0];
    for ( unsigned short i = 0; i < 2; ++i ) {
      gsl_spline2d_init( splines_[i], xa, ya, values_[i], q2_vals.size(), xbj_vals.size() );
    }
#else
    FatalError( Form( "GSL version ≥ 2.1 is required for bilinear interpolation.\n\tVersion %s is installed on this system!", GSL_VERSION ) );
#endif
  }

  GridHandler::~GridHandler()
  {
#ifdef GOOD_GSL
    for ( unsigned short i = 0; i < 2; ++i ) {
      gsl_spline2d_free( splines_[i] );
      if ( values_[i] ) delete[] values_[i];
    }
    if ( xacc_ ) gsl_interp_accel_free( xacc_ );
    if ( yacc_ ) gsl_interp_accel_free( yacc_ );
#endif
  }

  CepGen::StructureFunctions
  GridHandler::eval( double q2, double xbj ) const
  {
    CepGen::StructureFunctions ev;
#ifdef GOOD_GSL
    if ( gsl_spline2d_eval_e( splines_[0], q2, xbj, xacc_, yacc_, &ev.F2 ) == GSL_EDOM
      || gsl_spline2d_eval_e( splines_[1], q2, xbj, xacc_, yacc_, &ev.FL ) == GSL_EDOM ) {
      InWarning( Form( "Failed to evaluate the structure functions for Q² = %.5e GeV² / xbj = %.5e", q2, xbj ) );
      return ev;
    }
#else
    FatalError( Form( "GSL version ≥ 2.1 is required for bilinear interpolation.\n\tVersion %s is installed on this system!", GSL_VERSION ) );
#endif
    return ev;
  }
}
