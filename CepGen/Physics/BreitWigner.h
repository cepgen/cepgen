#ifndef CepGen_Physics_BreitWigner_h
#define CepGen_Physics_BreitWigner_h

#include <math.h>

namespace CepGen
{
  class BreitWigner
  {
    public:
      BreitWigner( double mean = 0., double gamma = 0., double emin = -1., double emax = -1. ) :
        mean_( mean ), gamma_( gamma ), emin_( emin ), emax_( emax ) {}
      BreitWigner( const BreitWigner& oth ) :
        mean_( oth.mean_ ), gamma_( oth.gamma_ ), emin_( oth.emin_ ), emax_( oth.emax_ ) {}

      double min() const { return emin_; }
      double max() const { return emax_; }
      inline double operator()( double x ) const {
        const double val = mean_+0.5*gamma_*tan( ( 2.*x-1. )*M_PI_2 );
        if ( ( emin_ >= 0. && val < emin_ ) || ( emax_ >= 0. && val > emax_ ) )
          return -1.;
        return val;
      }

    private:
      /// Mean of distribution
      double mean_;
      /// Width of distribution
      double gamma_;
      /// Minimal value
      double emin_;
      /// Maximal value
      double emax_;
  };
}

#endif
