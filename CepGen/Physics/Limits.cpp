#include "CepGen/Physics/Limits.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

namespace CepGen
{
  Limits::Limits( double min, double max ) :
    std::pair<double,double>( min, max )
  {}

  Limits::Limits( const Limits& rhs ) :
    std::pair<double,double>( rhs.first, rhs.second )
  {}

  void
  Limits::in( double low, double up )
  {
    first = low;
    second = up;
  }

  double
  Limits::range() const
  {
    if ( !hasMin() || !hasMax() )
      return 0.;
    return second-first;
  }

  bool
  Limits::hasMin() const
  {
    return first != kInvalid;
  }

  bool
  Limits::hasMax() const
  {
    return second != kInvalid;
  }

  bool
  Limits::passes( double val ) const
  {
    if ( hasMin() && val < min() )
      return false;
    if ( hasMax() && val > max() )
      return false;
    return true;
  }

  bool
  Limits::valid() const
  {
    return hasMin() || hasMax();
  }

  void
  Limits::save( bool& on, double& lmin, double& lmax ) const
  {
    on = false; lmin = lmax = 0.;
    if ( !valid() )
      return;
    on = true;
    if ( hasMin() )
      lmin = min();
    if ( hasMax() )
      lmax = max();
    if ( lmin == lmax )
      on = false;
  }

  double
  Limits::x( double v ) const
  {
    if ( v < 0. || v > 1. )
      CG_ERROR( "Limits:shoot" )
        << "x must be comprised between 0 and 1; x value = " << v << ".";
    if ( !valid() )
      return kInvalid;

    return first + ( second-first ) * v;
  }

  std::ostream&
  operator<<( std::ostream& os, const Limits& lim )
  {
    if ( !lim.hasMin() && !lim.hasMax() )
      return os << "no cuts";
    if ( !lim.hasMin() )
      return os << Form( "≤ %g", lim.max() );
    if ( !lim.hasMax() )
      return os << Form( "≥ %g", lim.min() );
    return os << Form( "%g → %g", lim.min(), lim.max() );
  }
}

