#include "CepGen/Core/GridParameters.h"
#include "CepGen/Core/Exception.h"
#include <cmath> // pow

namespace cepgen
{
  GridParameters::GridParameters( unsigned short ndim ) :
    gen_prepared( false ),
    f_max_diff( 0. ),
    r_boxes( 0 ),
    max_( pow( M_BIN, ndim ) ), f_max_( max_, 0. ), f_max_global_( 0. )
  {
    if ( ndim > MAX_DIM )
      throw CG_FATAL( "GridParameters" ) << "Phase space too large!\n\t"
        << "Either reduce the number of integration dimensions, or\n\t"
        << "increase the GridParameters::MAX_DIM parameter (not recommended).";

    //--- build and populate the grid
    coord_t coord( ndim );
    for ( unsigned int i = 0; i < max_; ++i ) {
      unsigned int jj = i;
      for ( unsigned int j = 0; j < ndim; ++j ) {
        unsigned int tmp = jj*INV_M_BIN;
        coord[j] = jj-tmp*M_BIN;
        jj = tmp;
      }
      n_map_.emplace_back( coord );
    }
    num.reserve( max_ );
  }

  const GridParameters::coord_t&
  GridParameters::n( size_t coord ) const
  {
    return n_map_.at( coord );
  }

  void
  GridParameters::setValue( size_t coord, double val )
  {
    //--- update function local and global maxima if needed
    f_max_[coord] = std::max( f_max_[coord], val );
    f_max_global_ = std::max( f_max_global_, val );
  }

  double
  GridParameters::maxValue( size_t coord ) const
  {
    return f_max_.at( coord );
  }

  void
  GridParameters::shoot( const gsl_rng* rng, size_t coord, std::vector<double>& out ) const
  {
    out.resize( 0 );
    for ( const auto& nv : n( coord ) )
      out.emplace_back( ( gsl_rng_uniform( rng )+nv ) * INV_M_BIN );
  }
}
