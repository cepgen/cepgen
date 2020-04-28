#include "CepGen/Utils/TimeKeeper.h"
#include "CepGen/Utils/String.h"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace cepgen
{
  namespace utils
  {
    TimeKeeper::TimeKeeper()
    {}

    void
    TimeKeeper::clear() {
      monitors_.clear();
      tmr_.reset();
    }

    TimeKeeper&
    TimeKeeper::tick( const std::string& func, double time )
    {
      monitors_[func].emplace_back( time > 0. ? time : tmr_.elapsed() );
      return *this;
    }

    std::string
    TimeKeeper::summary() const
    {
      if ( monitors_.empty() )
        return std::string();

      struct Monitor
      {
        std::string name;
        size_t size;
        double total, mean, rms;
        bool operator<( const Monitor& oth ) const { return total < oth.total; }
      };

      std::vector<Monitor> mons;
      for ( const auto& mon : monitors_ ) {
        const auto& tm = mon.second;
        const double total = tm.empty()
          ? -1.
          : std::accumulate( tm.begin(), tm.end(), 0. );
        const double mean = total/tm.size();
        const double rms = std::sqrt( std::fabs(
          std::inner_product( tm.begin(), tm.end(), tm.begin(), 0. )
            / tm.size() - mean*mean ) );
        mons.emplace_back( Monitor{ mon.first, tm.size(), total, mean, rms } );
      }
      std::sort( mons.rbegin(), mons.rend() );

      std::ostringstream oss;
      oss << utils::format( "%5s  %-100s   %12s\t%10s   %10s", "ncalls", "caller", "total (s)", "average (s)", "rms (s)" );
      for ( const auto& mon : mons )
        oss
          << utils::format( "\n* [%10u | %-100s] %12.6f\t%10e +/- %10e", mon.size, mon.name.c_str(), mon.total, mon.mean, mon.rms );

      return oss.str();
    }

    TimeKeeper::Ticker::Ticker( TimeKeeper& tk, const std::string& name ) :
      tk_( tk ), name_( name )
    {}

    TimeKeeper::Ticker::~Ticker()
    {
      tk_.tick( name_, tmr_.elapsed() );
    }
  }
}
