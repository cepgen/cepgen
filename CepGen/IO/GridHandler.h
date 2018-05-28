#ifndef CepGen_IO_GridHandler_h
#define CepGen_IO_GridHandler_h

#include <gsl/gsl_version.h>
#ifdef GSL_MAJOR_VERSION
#  if GSL_MAJOR_VERSION > 2 || ( GSL_MAJOR_VERSION == 2 && GSL_MINOR_VERSION >= 1 )
#    include <gsl/gsl_interp2d.h>
#    include <gsl/gsl_spline2d.h>
#    define GOOD_GSL 1
#  endif
#endif
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include "CepGen/Core/Exception.h"
#include <set>

namespace CepGen
{
  class StructureFunctions;
  enum struct GridType
  {
    kLinear = 0,
    kLogarithmic = 1,
    kSquare = 2
  };
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
      explicit GridHandler( const GridType& grid_type ) :
#ifdef GOOD_GSL
        xacc_( gsl_interp_accel_alloc(), gsl_interp_accel_free ),
        yacc_( gsl_interp_accel_alloc(), gsl_interp_accel_free ),
#endif
        grid_type_( grid_type )
      {}
      ~GridHandler() {
#ifdef GOOD_GSL
        for ( unsigned short i = 0; i < splines_.size(); ++i ) {
          gsl_spline2d_free( splines_[i] );
          if ( values_[i] )
            delete[] values_[i];
        }
#endif
      }

      /// Interpolate a point to a given coordinate
      grid_t eval( double x, double y ) const {
        grid_t ev;
        if ( grid_type_ == GridType::kLogarithmic ) {
          x = log10( x );
          y = log10( y );
        }
        else if ( grid_type_ == GridType::kSquare ) {
          x *= x;
          y *= y;
        }
#ifdef GOOD_GSL
        unsigned short i = 0;
        for ( auto& sp : splines_ )
          if ( gsl_spline2d_eval_e( sp, x, y, xacc_.get(), yacc_.get(), &ev.value[i++] ) != GSL_SUCCESS )
            CG_WARNING( "GridHandler" )
              << "Failed to evaluate the grid value "
              << "for x = " << x << " / y = " << y << ".";
#else
        double x_1 = -1., x_2 = -1., y_1 = -1., y_2 = -1.;
        std::array<double,N> ext_11, ext_12, ext_21, ext_22;
        for ( unsigned short i = 0; i < x_vals_.size()-1; ++i ) {
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
          if ( x_1 == val.x && y_1 == val.y )
            for ( unsigned short i = 0; i < ext_11.size(); ++i )
              ext_11[i] = val.value[i];
          if ( x_2 == val.x && y_1 == val.y )
            for ( unsigned short i = 0; i < ext_12.size(); ++i )
              ext_12[i] = val.value[i];
          if ( x_1 == val.x && y_2 == val.y )
            for ( unsigned short i = 0; i < ext_21.size(); ++i )
              ext_21[i] = val.value[i];
          if ( x_2 == val.x && y_2 == val.y )
            for ( unsigned short i = 0; i < ext_22.size(); ++i )
              ext_22[i] = val.value[i];
        }
        // now that we have the boundaries, we may interpolate
        const double x2x1 = x_2-x_1;
        const double y2y1 = y_2-y_1;
        const double x2x = x_2-x;
        const double y2y = y_2-y;
        const double xx1 = x-x_1;
        const double yy1 = y-y_1;
        for ( unsigned short i = 0; i < ext_11.size(); ++i )
          ev.value[i] = 1.0 / ( x2x1*y2y1 ) * ( ext_11[i]*x2x*y2y + ext_21[i]*xx1*y2y + ext_12[i]*x2x*yy1 + ext_22[i]*xx1*yy1 );
#endif
        return ev;
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
      std::unique_ptr<gsl_interp_accel,void(*)( gsl_interp_accel* )> xacc_, yacc_;
      std::array<double*,N> values_;
#else
      /// List of coordinates handled in the grid
      std::vector<double> x_vals_;
      /// List of coordinates handled in the grid
      std::vector<double> y_vals_;
#endif
      GridType grid_type_;

  };
}

#endif

