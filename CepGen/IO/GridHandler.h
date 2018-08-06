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
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include "CepGen/Core/Exception.h"
#include <vector>

namespace CepGen
{
  class StructureFunctions;
  enum struct GridType
  {
    linear = 0,
    logarithmic = 1,
    square = 2
  };
  /// A generic class for D-dimensional grid interpolation
  /// \param N Number of values handled per point
  template <size_t D,size_t N=1>
  class GridHandler
  {
    public:
      typedef std::array<double,D> coord_t;
      typedef std::array<double,N> values_t;
      /// A point and its associated value(s)
      struct point_t
      {
        /// Coordinate(s)
        coord_t coord;
        /// (list of) value(s) for this point
        values_t value;
      };

    public:
      explicit GridHandler( const GridType& grid_type ) :
        grid_type_( grid_type ), accel_{}
      {
        for ( size_t i = 0; i < D; ++i )
          accel_.emplace_back( std::unique_ptr<gsl_interp_accel,void(*)( gsl_interp_accel* )>( gsl_interp_accel_alloc(), gsl_interp_accel_free ) );
      }
      ~GridHandler() {
#ifdef GOOD_GSL
        for ( size_t i = 0; i < D; ++i )
          if ( values_[i] )
            delete[] values_[i];
#endif
      }

      /// Interpolate a point to a given coordinate
      values_t eval( coord_t in_coords ) const {
        values_t out;
        coord_t coord = in_coords;
        switch ( grid_type_ ) {
          case GridType::logarithmic: {
            for ( auto& c : coord )
              c = log10( c );
          } break;
          case GridType::square: {
            for ( auto& c : coord )
              c *= c;
          } break;
          default: break;
        }
        //--- dimension of the vector space coordinate to evaluate
        switch ( D ) {
          case 1: {
            for ( size_t i = 0; i < N; ++i ) {
              int res = gsl_spline_eval_e( splines_1d_.at( i ).get(), coord.at( 0 ), accel_.at( 0 ).get(), &out[i] );
              if ( res != GSL_SUCCESS ) {
                out[i] = 0.;
                CG_WARNING( "GridHandler" )
                  << "Failed to evaluate the grid value (N=" << i << ") "
                  << "for x = " << in_coords.at( 0 ) << ". "
                  << "GSL error: " << gsl_strerror( res );
              }
            }
          } break;
          case 2: {
            const double x = coord.at( 0 ), y = coord.at( 1 );
#ifdef GOOD_GSL
            for ( size_t i = 0; i < N; ++i ) {
              int res = gsl_spline2d_eval_e( splines_2d_.at( i ).get(), x, y, accel_.at( 0 ).get(), accel_.at( 1 ).get(), &out[i] );
              if ( res != GSL_SUCCESS ) {
                out[i] = 0.;
                CG_WARNING( "GridHandler" )
                  << "Failed to evaluate the grid value (N=" << i << ") "
                  << "for x = " << in_coords.at( 0 ) << " / y = " << in_coords.at( 1 ) << ". "
                  << "GSL error: " << gsl_strerror( res );
              }
            }
#else
            //--- retrieve the indices of the bin in the set
            coord_t min, max;
            findIndices( coord, min, max );
            //--- find boundaries values
            gridpoint_t ext_11, ext_12, ext_21, ext_22;
            for ( const auto& val : values_raw_ ) {
              if      ( min[0] == val.coord.at( 0 ) && min[1] == val.coord.at( 1 ) ) ext_11 = val.value;
              else if ( min[0] == val.coord.at( 0 ) && max[1] == val.coord.at( 1 ) ) ext_12 = val.value;
              else if ( max[0] == val.coord.at( 0 ) && min[1] == val.coord.at( 1 ) ) ext_21 = val.value;
              else if ( max[0] == val.coord.at( 0 ) && max[1] == val.coord.at( 1 ) ) ext_22 = val.value;
            }
            //--- now that we have the boundaries, we may interpolate
            coord_t c_d;
            for ( size_t i = 0; i < D; ++i )
              c_d[i] = ( coord.at( i )-min[i] )/( max[i]-min[i] );
            const gridpoint_t ext_1 = ext_11*( 1.-c_d[0] ) + ext_21*c_d[0];
            const gridpoint_t ext_2 = ext_12*( 1.-c_d[0] ) + ext_22*c_d[0];
            out = ext_1*( 1.-c_d[1] )+ext_2*c_d[1];
#endif
          } break;
          case 3: {
            //--- retrieve the indices of the bin in the set
            coord_t min, max;
            findIndices( coord, min, max );
            //--- find boundaries values
            gridpoint_t ext_111, ext_112, ext_121, ext_122, ext_211, ext_212, ext_221, ext_222;
            for ( const auto& val : values_raw_ ) {
              if      ( min[0] == val.coord.at( 0 ) && min[1] == val.coord.at( 1 ) && min[2] == val.coord.at( 2 ) ) ext_111 = val.value;
              else if ( min[0] == val.coord.at( 0 ) && min[1] == val.coord.at( 1 ) && max[2] == val.coord.at( 2 ) ) ext_112 = val.value;
              else if ( min[0] == val.coord.at( 0 ) && max[1] == val.coord.at( 1 ) && min[2] == val.coord.at( 2 ) ) ext_121 = val.value;
              else if ( min[0] == val.coord.at( 0 ) && max[1] == val.coord.at( 1 ) && max[2] == val.coord.at( 2 ) ) ext_122 = val.value;
              else if ( max[0] == val.coord.at( 0 ) && min[1] == val.coord.at( 1 ) && min[2] == val.coord.at( 2 ) ) ext_211 = val.value;
              else if ( max[0] == val.coord.at( 0 ) && min[1] == val.coord.at( 1 ) && max[2] == val.coord.at( 2 ) ) ext_212 = val.value;
              else if ( max[0] == val.coord.at( 0 ) && max[1] == val.coord.at( 1 ) && min[2] == val.coord.at( 2 ) ) ext_221 = val.value;
              else if ( max[0] == val.coord.at( 0 ) && max[1] == val.coord.at( 1 ) && max[2] == val.coord.at( 2 ) ) ext_222 = val.value;
            }
            //--- now that we have the boundaries, we may interpolate
            coord_t c_d;
            for ( size_t i = 0; i < D; ++i )
              c_d[i] = ( coord.at( i )-min[i] )/( max[i]-min[i] );
            const gridpoint_t ext_11 = ext_111*( 1.-c_d[0] ) + ext_211*c_d[0];
            const gridpoint_t ext_12 = ext_112*( 1.-c_d[0] ) + ext_212*c_d[0];
            const gridpoint_t ext_21 = ext_121*( 1.-c_d[0] ) + ext_221*c_d[0];
            const gridpoint_t ext_22 = ext_122*( 1.-c_d[0] ) + ext_222*c_d[0];
            const gridpoint_t ext_1 = ext_11*( 1.-c_d[1] ) + ext_21*c_d[1];
            const gridpoint_t ext_2 = ext_12*( 1.-c_d[1] ) + ext_22*c_d[1];
            out = ext_1*( 1.-c_d[2] )+ext_2*c_d[2];
          } break;
          default:
            throw CG_FATAL( "GridHandler" ) << "Unsupported number of dimensions: " << N;
        }
        return out;
      }

      /// Insert a new value in the grid
      void insert( const point_t& point ) {
        auto mod_coord = point.coord;
        if ( grid_type_ != GridType::linear )
          for ( auto& c : mod_coord )
            switch ( grid_type_ ) {
              case GridType::logarithmic:
                c = log10( c ); break;
              case GridType::square:
                c *= c; break;
              default: break;
            }
        values_raw_.emplace_back( point_t{ mod_coord, point.value } );
      }
      /// Return the list of values handled in the grid
      std::vector<point_t> values() const { return values_raw_; }

      /// Initialise the grid and all useful interpolators/accelerators
      void init() {
        if ( values_raw_.empty() )
          CG_ERROR( "GridHandler" ) << "Empty grid.";
        gsl_set_error_handler_off();
        //--- start by building grid coordinates from raw values
        for ( auto& c : coords_ )
          c.clear();
        for ( const auto& val : values_raw_ ) {
          unsigned short i = 0;
          for ( const auto& c : val.coord ) {
            if ( std::find( coords_.at( i ).begin(), coords_.at( i ).end(), c ) == coords_.at( i ).end() )
              coords_.at( i ).emplace_back( c );
            ++i;
          }
        }
        for ( auto& c : coords_ )
          std::sort( c.begin(), c.end() );
        { //--- debugging of the grid coordinates
          std::ostringstream os;
          unsigned short i = 0;
          for ( const auto& cs : coords_ ) {
            os << "coordinate " << (i++) << " has " << cs.size() << " member(s):";
            unsigned short j = 0;
            for ( const auto& val : cs )
              os << ( j++ % 20 == 0 ? "\n  " : " " ) << val;
            os << "\n";
          }
          CG_DEBUG( "GridHandler" ) << "Grid dump:\n" << os.str();
        }
        //--- particularise by dimension
        switch ( D ) {
          case 1: { //--- x |-> (f1,...)
            const gsl_interp_type* type = gsl_interp_cspline;
            //const gsl_interp_type* type = gsl_interp_steffen;
            if ( gsl_interp_type_min_size( type ) >= values_raw_.size() )
              throw CG_FATAL( "GridHandler" ) << "Not enough points for \"" << type->name << "\" type of interpolation.\n\t"
                << "Minimum required: " << gsl_interp_type_min_size( type ) << ", got " << values_raw_.size() << "!";
            for ( size_t i = 0; i < N; ++i ) {
              values_[i] = new double[values_raw_.size()];
              splines_1d_.emplace_back( gsl_spline_alloc( type, values_raw_.size() ), gsl_spline_free );
            }
            std::vector<double> x_vec;
            for ( unsigned short i = 0; i < values_raw_.size(); ++i ) {
              const auto& val = values_raw_.at( i );
              x_vec.emplace_back( val.coord.at( 0 ) );
              for ( size_t j = 0; j < N; ++j )
                values_[j][i] = val.value.at( j );
            }
            for ( unsigned short i = 0; i < splines_1d_.size(); ++i )
              gsl_spline_init( splines_1d_.at( i ).get(), &x_vec[0], values_.at( i ), values_raw_.size() );
          } break;
          case 2: { //--- (x,y) |-> (f1,...)
#ifdef GOOD_GSL
            const gsl_interp2d_type* type = gsl_interp2d_bilinear;
            splines_2d_.clear();
            for ( size_t i = 0; i < N; ++i ) {
              values_[i] = new double[coords_.at( 0 ).size() * coords_.at( 1 ).size()];
              splines_2d_.emplace_back( gsl_spline2d_alloc( type, coords_.at( 0 ).size(), coords_.at( 1 ).size() ), gsl_spline2d_free );
            }

            // second loop over all points to populate the grid
            for ( const auto& val : values_raw_ ) {
              double val_x = val.coord.at( 0 ), val_y = val.coord.at( 1 );
              // retrieve the index of the bin in the set
              const unsigned short id_x = std::distance( coords_.at( 0 ).begin(), std::lower_bound( coords_.at( 0 ).begin(), coords_.at( 0 ).end(), val_x ) );
              const unsigned short id_y = std::distance( coords_.at( 1 ).begin(), std::lower_bound( coords_.at( 1 ).begin(), coords_.at( 1 ).end(), val_y ) );
              for ( unsigned short i = 0; i < splines_2d_.size(); ++i )
                gsl_spline2d_set( splines_2d_.at( i ).get(), values_[i], id_x, id_y, val.value[i] );
            }

            // initialise splines objects
            std::vector<double> x_vec( coords_.at( 0 ).begin(), coords_.at( 0 ).end() ), y_vec( coords_.at( 1 ).begin(), coords_.at( 1 ).end() );
            for ( unsigned short i = 0; i < splines_2d_.size(); ++i )
              gsl_spline2d_init( splines_2d_.at( i ).get(), &x_vec[0], &y_vec[0], values_.at( i ), x_vec.size(), y_vec.size() );
#else
            CG_WARNING( "GridHandler" )
              << "GSL version ≥ 2.1 is required for spline bilinear interpolation.\n\t"
              << "Version " << GSL_VERSION << " is installed on this system!\n\t"
              << "Will use a linear approximation instead.";
#endif
          } break;
        }
      }

    protected:
      GridType grid_type_;
      /// List of coordinates and associated value(s) in the grid
      std::vector<point_t> values_raw_;

      std::vector<std::unique_ptr<gsl_interp_accel,void(*)( gsl_interp_accel* )> > accel_;
      std::vector<std::unique_ptr<gsl_spline,void(*)( gsl_spline* )> > splines_1d_;
#ifdef GOOD_GSL
      std::vector<std::unique_ptr<gsl_spline2d,void(*)( gsl_spline2d* )> > splines_2d_;
#endif
      std::array<std::vector<double>,D> coords_;
      std::array<double*,N> values_;

    private:
      void findIndices( const coord_t& coord, coord_t& min, coord_t& max ) const {
        for ( size_t i = 0; i < D; ++i ) {
          if ( coord.at( i ) < *coords_.at( i ).begin() )
            min[i] = max[i] = *coords_.at( i ).begin();
          else if ( coord.at( i ) > *coords_.at( i ).rbegin() )
            min[i] = max[i] = *coords_.at( i ).rbegin();
          else {
            auto it_coord = std::lower_bound( coords_.at( i ).begin(), coords_.at( i ).end(), coord.at( i ) );
            const unsigned short id = std::distance( coords_.at( i ).begin(), it_coord );
            min[i] = coords_.at( i ).at( id );
            max[i] = coords_.at( i ).at( id < coords_.at( i ).size()-1 ? id+1 : id );
          }
        }
      }

      struct gridpoint_t : values_t
      {
        gridpoint_t() : values_t() {}
        gridpoint_t( const values_t& arr ) : values_t( arr ) {}
        gridpoint_t operator*( double c ) const {
          gridpoint_t out = *this;
          for ( auto& a : out )
            a *= c;
          return out;
        }
        gridpoint_t operator+( const gridpoint_t& rhs ) const {
          gridpoint_t out = *this;
          for ( size_t i = 0; i < out.size(); ++i )
            out[i] += rhs[i];
          return out;
        }
      };
  };
}

#endif

