#include "CepGen/Utils/ProgressBar.h"

#include <iomanip>
#include <iterator>
#include <algorithm>
#include <sstream>

#include <math.h>

namespace cepgen
{
  namespace utils
  {
    ProgressBar::ProgressBar( size_t tot, size_t freq ) :
      bar_pattern_( BAR_LENGTH, '=' ), total_( tot ), frequency_( freq )
    {}

    void ProgressBar::update( size_t iter ) const
    {
      const size_t percent = iter*100./total_;
      if ( percent % frequency_ == 0 || iter == total_ ) {
        int lpad = int( percent/100. * BAR_LENGTH );
        int rpad = BAR_LENGTH-lpad;
        fprintf( stderr, "\r%3zu%% [%.*s%*s]", percent, lpad, bar_pattern_.c_str(), rpad, "" );
        fflush( stderr );
        if ( iter == total_ )
          fprintf( stderr, "\n" );
      }
    }
  }
}
