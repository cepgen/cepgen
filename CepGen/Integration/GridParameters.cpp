#include "CepGen/Integration/GridParameters.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Core/Exception.h"

#include <cmath> // pow

namespace cepgen
{
  GridParameters::GridParameters( size_t ndim ) :
    gen_prepared( false ),
    f_max_diff( 0. ),
    ndim_( ndim ), f_max_global_( 0. )
  {
    //--- build and populate the grid
    coord_t coord( ndim, 0 );
    for ( size_t i = 0; i < pow( M_BIN, ndim_ ); ++i ) {
      size_t jj = i;
      for ( size_t j = 0; j < ndim; ++j ) {
        size_t tmp = jj*INV_M_BIN;
        coord[j] = jj-tmp*M_BIN;
        jj = tmp;
      }
      n_map_.emplace_back( point_t{ coord, 0ul, 0. } );
    }
  }

  size_t
  GridParameters::size() const
  {
    return n_map_.size();
  }

  const GridParameters::coord_t&
  GridParameters::n( size_t coord ) const
  {
    return n_map_.at( coord ).coordinates;
  }

  void
  GridParameters::setValue( size_t coord, double val )
  {
    //--- update function local and global maxima if needed
    n_map_.at( coord ).f_max = std::max( n_map_.at( coord ).f_max, val );
    f_max_global_ = std::max( f_max_global_, val );
  }

  double
  GridParameters::maxValue( size_t coord ) const
  {
    return n_map_.at( coord ).f_max;
  }

  size_t
  GridParameters::numPoints( size_t coord ) const
  {
    return n_map_.at( coord ).num_points;
  }

  void
  GridParameters::increment( size_t coord )
  {
    n_map_.at( coord ).num_points++;
  }

  void
  GridParameters::shoot( const Integrator* integr, size_t coord, std::vector<double>& out ) const
  {
    const auto& nv = n_map_.at( coord );
    for ( size_t i = 0; i < nv.coordinates.size(); ++i )
      out[i] = ( integr->uniform()+nv.coordinates.at( i ) ) * INV_M_BIN;
  }

  void
  GridParameters::dump() const
  {
    CG_INFO( "GridParameters:dump" ).log( [&]( auto& info ) {
      size_t i = 0;
      for ( const auto& point : n_map_ )
        info << "\nn[" << i++ << "]: "
          << "coord=" << point.coordinates << ", "
          << "num points: " << point.num_points << ", "
          << "max=" << point.f_max << ".";
    } );
  }
}
