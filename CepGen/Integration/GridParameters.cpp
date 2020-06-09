#include "CepGen/Integration/GridParameters.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Core/Exception.h"

#include <cmath> // pow

namespace cepgen
{
  GridParameters::GridParameters( size_t ndim ) :
    gen_prepared( false ),
    f_max_diff( 0. ),
    max_( pow( M_BIN, ndim ) ), num_points_( max_, 0ul ),
    f_max_( max_, 0. ), f_max_global_( 0. )
  {
    if ( ndim > MAX_DIM )
      throw CG_FATAL( "GridParameters" ) << "Phase space too large!\n\t"
        << "Either reduce the number of integration dimensions, or\n\t"
        << "increase the GridParameters::MAX_DIM parameter (not recommended).";

    //--- build and populate the grid
    coord_t coord( ndim );
    for ( size_t i = 0; i < max_; ++i ) {
      size_t jj = i;
      for ( size_t j = 0; j < ndim; ++j ) {
        size_t tmp = jj*INV_M_BIN;
        //coord[j] = roundf( jj-tmp*M_BIN );
        coord[j] = jj-tmp*M_BIN;
        jj = tmp;
      }
      n_map_.emplace_back( coord );
    }
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
    f_max_[coord] = std::max( f_max_.at( coord ), val );
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
  GridParameters::setTrial( size_t coord )
  {
    num_points_[coord]++;
  }

  void
  GridParameters::shoot( const Integrator* integr, size_t coord, std::vector<double>& out ) const
  {
    const auto& nv = n_map_.at( coord );
    for ( size_t i = 0; i < nv.size(); ++i )
      out[i] = ( integr->uniform()+nv.at( i ) ) * INV_M_BIN;
  }

  void
  GridParameters::dump() const
  {
    std::ostringstream os;
    size_t i = 0;
    for ( const auto& n : n_map_ ) {
      os << "n[" << i++ << "] = {";
      std::string sep;
      for ( const auto& v : n )
        os << sep << v, sep = ", ";
      os << "}" << std::endl;
    }
    CG_INFO( "GridParameters:dump" ) << os.str();
  }
}
