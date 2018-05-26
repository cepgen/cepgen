#ifndef CepGen_IO_GridHandler_h
#define CepGen_IO_GridHandler_h

#include <gsl/gsl_version.h>

#ifdef GSL_MAJOR_VERSION
#if GSL_MAJOR_VERSION > 2 || ( GSL_MAJOR_VERSION == 2 && GSL_MINOR_VERSION >= 1 )
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#define GOOD_GSL 1
#endif
#endif

#include <array>
#include <vector>
#include <set>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <fstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

namespace CepGen
{
  class StructureFunctions;
  /// A generic class for 2-dimensional grid interpolation
  /// \param N Number of values handled per point
  template <size_t N=1>
  class GridHandler
  {
    public:
      /// A point and its associated value(s)
      struct grid_t
      {
        /// First coordinate
        float x;
        /// Second coordinate
        float y;
        /// (list of) value(s) for this point
        double value[N];
      };

    public:
      explicit GridHandler()
#ifdef GOOD_GSL
        : xacc_( nullptr ), yacc_( nullptr )
#endif
      {}
      ~GridHandler() {
#ifdef GOOD_GSL
        for ( unsigned short i = 0; i < splines_.size(); ++i ) {
          gsl_spline2d_free( splines_[i] );
          if ( values_[i] )
            delete[] values_[i];
        }
        if ( xacc_ )
          gsl_interp_accel_free( xacc_ );
        if ( yacc_ )
          gsl_interp_accel_free( yacc_ );
#endif
      }

      /// Interpolate a point to a given coordinate
      grid_t eval( double x, double y ) const {
        grid_t ev;
#ifdef GOOD_GSL
        unsigned short i = 0;
        for ( auto& sp : splines_ )
          if ( gsl_spline2d_eval_e( sp, log10( x ), log10( y ), xacc_, yacc_, &ev.value[i++] ) != GSL_SUCCESS )
            CG_WARNING( "GridHandler" )
              << "Failed to evaluate the grid value "
              << "for x = " << x << " / y = " << y << ".";
        return ev;
#else
        double x_1 = -1., x_2 = -1., y_1 = -1., y_2 = -1.;
        double f2_11 = 0., f2_12 = 0., f2_21 = 0., f2_22 = 0.;
        double fl_11 = 0., fl_12 = 0., fl_21 = 0., fl_22 = 0.;
        for ( unsigned int i = 0; i < x_vals_.size()-1; ++i ) {
          if ( x >= x_vals_.at( i ) && x < x_vals_.at( i+1 ) ) {
            x_1 = x_vals_.at( i );
            x_2 = x_vals_.at( i+1 );
            break;
          }
        }
        for ( unsigned int i = 0; i < y_vals_.size()-1; ++i ) {
          if ( y >= y_vals_.at( i ) && y < y_vals_.at( i+1 ) ) {
            y_1 = y_vals_.at( i );
            y_2 = y_vals_.at( i+1 );
            break;
          }
        }
        for ( const auto& val : values_raw_ ) {
          if ( x_1 == val.x && y_1 == val.y ) {
            f2_11 = val.f2;
            fl_11 = val.fl;
          }
          if ( x_2 == val.x && y_1 == val.y ) {
            f2_12 = val.f2;
            fl_12 = val.fl;
          }
          if ( x_1 == val.x && y_2 == val.y ) {
            f2_21 = val.f2;
            fl_21 = val.fl;
          }
          if ( x_2 == val.x && y_2 == val.y ) {
            f2_22 = val.f2;
            fl_22 = val.fl;
          }
        }
        const double x2x1 = y_2-y_1;
        const double y2y1 = x_2-x_1;
        const double x2x = y_2-y;
        const double y2y = x_2-x;
        const double yy1 = x-x_1;
        const double xx1 = y-y_1;
        ev.F2 = 1.0 / ( x2x1*y2y1 ) * ( f2_11*x2x*y2y + f2_21*xx1*y2y + f2_12*x2x*yy1 + f2_22*xx1*yy1 );
        ev.FL = 1.0 / ( x2x1*y2y1 ) * ( fl_11*x2x*y2y + fl_21*xx1*y2y + fl_12*x2x*yy1 + fl_22*xx1*yy1 );
        return ev;
#endif
      }

      /// Return the list of values handled in the grid
      std::vector<grid_t> values() const { return values_raw_; }

    protected:
      /// Initialise the grid and all useful interpolators/accelerators
      void initGSL( const std::set<double>& x_vals, const std::set<double>& y_vals ) {
#ifdef GOOD_GSL
        gsl_set_error_handler_off();
        const gsl_interp2d_type* type = gsl_interp2d_bilinear;

        for ( unsigned short i = 0; i < splines_.size(); ++i ) {
          values_[i] = new double[x_vals.size() * y_vals.size()];
          splines_[i] = gsl_spline2d_alloc( type, x_vals.size(), y_vals.size() );
        }
        xacc_ = gsl_interp_accel_alloc();
        yacc_ = gsl_interp_accel_alloc();

        // second loop over all points to populate the grid
        for ( const auto& val : values_raw_ ) {
          // retrieve the index of the Q2/xbj bin in the set
          const unsigned short id_x = std::distance( x_vals.begin(), x_vals.lower_bound( log10( val.x  ) ) ),
                               id_y = std::distance( y_vals.begin(), y_vals.lower_bound( log10( val.y ) ) );
          unsigned short i = 0;
          for ( auto& sp : splines_ ) {
            gsl_spline2d_set( sp, values_[i], id_x, id_y, val.value[i] );
            ++i;
          }
        }

        // initialise splines objects
        std::vector<double> q2_vec( x_vals.begin(), x_vals.end() ), xbj_vec( y_vals.begin(), y_vals.end() );
        for ( unsigned short i = 0; i < splines_.size(); ++i )
          gsl_spline2d_init( splines_[i], &q2_vec[0], &xbj_vec[0], values_[i], x_vals.size(), y_vals.size() );
#else
        CG_WARNING( "GridHandler" )
          << "GSL version ≥ 2.1 is required for spline bilinear interpolation.\n\t"
          << "Version " << GSL_VERSION << " is installed on this system!\n\t"
          << "Will use a linear approximation instead.\n\t"
          << "You may check the numerical validity of this approach...";
#endif
      }

      /// List of coordinates and associated value(s) in the grid
      std::vector<grid_t> values_raw_;
#ifdef GOOD_GSL
      std::array<gsl_spline2d*,N> splines_;
      gsl_interp_accel* xacc_, *yacc_;
      std::array<double*,N> values_;
#else
      /// List of coordinates handled in the grid
      std::vector<double> x_vals_;
      /// List of coordinates handled in the grid
      std::vector<double> y_vals_;
#endif

  };
}

#endif

