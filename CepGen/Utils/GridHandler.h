#ifndef CepGen_Core_GridHandler_h
#define CepGen_Core_GridHandler_h

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

#include "CepGen/Core/Exception.h"
#include <vector>
#include <map>

namespace cepgen
{
  namespace strfun { class Parameterisation; }
  /// Interpolation type for the grid coordinates
  enum struct GridType { linear, logarithmic, square };
  /// \brief A generic class for D-dimensional grid interpolation
  /// \tparam D Number of variables in the grid (dimension)
  /// \tparam N Number of values handled per point
  template <size_t D,size_t N=1>
  class GridHandler
  {
    public:
      typedef std::vector<double> coord_t; ///< Coordinates container
      typedef std::array<double,N> values_t; ///< Value(s) at a given coordinate

    public:
      /// Build a grid interpolator from a grid type
      explicit GridHandler( const GridType& grid_type );
      ~GridHandler() {}

      /// Interpolate a point to a given coordinate
      values_t eval( coord_t in_coords ) const;

      /// Insert a new value in the grid
      void insert( coord_t coord, values_t value );
      /// Return the list of values handled in the grid
      std::map<coord_t,values_t> values() const { return values_raw_; }

      /// Initialise the grid and all useful interpolators/accelerators
      void init();
      /// Grid boundaries (collection of pair(min,max))
      std::array<std::pair<double,double>,D> boundaries() const;

    protected:
      /// Type of interpolation for the grid members
      GridType grid_type_;
      /// List of coordinates and associated value(s) in the grid
      std::map<coord_t,values_t> values_raw_;
      /// GSL grid interpolation accelerator
      std::vector<std::unique_ptr<gsl_interp_accel,void(*)( gsl_interp_accel* )> > accel_;
      /// Collection of splines for linear interpolations
      std::vector<std::unique_ptr<gsl_spline,void(*)( gsl_spline* )> > splines_1d_;
#ifdef GOOD_GSL
      /// Collection of splines for bilinear interpolations
      std::vector<std::unique_ptr<gsl_spline2d,void(*)( gsl_spline2d* )> > splines_2d_;
#endif
      /// Collection of coordinates building up the grid
      std::array<coord_t,D> coords_;
      /// Collection of values for all points in the grid
      std::array<std::unique_ptr<double[]>,N> values_;

    private:
      /// Retrieve lower and upper grid indices for a given coordinate
      void findIndices( const coord_t& coord, coord_t& min, coord_t& max ) const;
      /// A single value in grid coordinates
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
  template class GridHandler<3,1>;
  template class GridHandler<2,2>;
}

#endif
