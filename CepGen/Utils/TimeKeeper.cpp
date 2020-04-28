#include "CepGen/Utils/TimeKeeper.h"

#include <numeric>
#include <sstream>
#include <cmath>

namespace cepgen
{
  namespace utils
  {
    TimeKeeper::TimeKeeper()
    {}

    void
    TimeKeeper::clear() {
      timers_.clear();
      tmr_.reset();
    }

    TimeKeeper&
    TimeKeeper::tick( const std::string& func )
    {
      timers_[func].emplace_back( tmr_.elapsed() );
      return *this;
    }

    std::string
    TimeKeeper::summary() const
    {
      std::ostringstream oss;
      std::string sep;
      for ( const auto& tmr : timers_ ) {
        const double mean = tmr.second.empty()
          ? -1.
          : std::accumulate( tmr.second.begin(), tmr.second.end(), 0. )/tmr.second.size();
        const double rms = std::sqrt(
          std::inner_product( tmr.second.begin(), tmr.second.end(), tmr.second.begin(), 0. )
            /tmr.second.size()-mean*mean );
        oss << sep
          << "[" << tmr.first << ":" << tmr.second.size() << "] "
          << mean << " +/- " << rms, sep = "\n";
      }
      return oss.str();
    }
  }
}
