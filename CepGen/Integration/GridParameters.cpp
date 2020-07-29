#include "CepGen/Integration/GridParameters.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Core/Exception.h"

#include <cmath> // pow

namespace cepgen
{
  GridParameters::GridParameters( size_t ndim ) :
    gen_prepared( false ),
    correc( 0. ), correc2( 0. ),
    f_max2( 0. ), f_max_diff( 0. ), f_max_old( 0. ),
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
      coords_.emplace_back( coord );
      num_points_.emplace_back( 0ul );
      f_max_.emplace_back( 0. );
    }
  }

  size_t
  GridParameters::size() const
  {
    return coords_.size();
  }

  const GridParameters::coord_t&
  GridParameters::n( size_t coord ) const
  {
    return coords_.at( coord );
  }

  void
  GridParameters::setValue( size_t coord, double val )
  {
    //--- update function local and global maxima if needed
    f_max_.at( coord ) = std::max( f_max_.at( coord ), val );
    f_max_global_ = std::max( f_max_global_, val );
  }

  double
  GridParameters::maxValue( size_t coord ) const
  {
    return f_max_.at( coord );
  }

  size_t
  GridParameters::numPoints( size_t coord ) const
  {
    return num_points_.at( coord );
  }

  void
  GridParameters::increment( size_t coord )
  {
    num_points_.at( coord )++;
  }

  void
  GridParameters::shoot( const Integrator* integr, size_t coord, std::vector<double>& out ) const
  {
    const auto& nv = coords_.at( coord );
    for ( size_t i = 0; i < nv.size(); ++i )
      out[i] = ( integr->uniform()+nv.at( i ) ) * INV_M_BIN;
  }

  void
  GridParameters::dump() const
  {
    CG_INFO( "GridParameters:dump" ).log( [&]( auto& info ) {
      for ( size_t i = 0; i < coords_.size(); ++i )
        info << "\nn[" << i << "]: "
          << "coord=" << coords_.at( i ) << ", "
          << "num points: " << num_points_.at( i ) << ", "
          << "max=" << f_max_.at( i ) << ".";
    } );
  }
}
