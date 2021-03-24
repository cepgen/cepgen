#include "CepGen/Utils/TimeKeeper.h"
#include "CepGen/Utils/String.h"

#include <algorithm>
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

      //--- a bit of arithmetics

      struct Monitor
      {
        std::string name;
        size_t size;
        double total, mean, rms;
        bool operator<( const Monitor& oth ) const { return total < oth.total; }
      };

      std::vector<Monitor> mons;
      double total_time = 0.;
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
        total_time += total;
      }

      //--- sort by total clock time (desc.)

      std::sort( mons.rbegin(), mons.rend() );

      //--- displaying the various probes

      std::ostringstream oss;
      oss << utils::format( "%2s | %-90s | %12s\t%10s\t%5s",
                            "#", "Caller", "Total (ms)",
                            "Average (ms)", "RMS (ms)" );
      for ( const auto& mon : mons )
        oss
          << utils::format( "\n%10u | %-90s | %12.6f\t%10e\t%5.3e",
                            mon.size, mon.name.c_str(), mon.total*1.e3,
                            mon.mean*1.e3, mon.rms*1.e3 );

      return oss.str();
    }

    TimeKeeper::Ticker::Ticker( TimeKeeper* tk, const std::string& name ) :
      tk_( tk ), name_( name )
    {}

    TimeKeeper::Ticker::~Ticker()
    {
      if ( tk_ )
        tk_->tick( name_, tmr_.elapsed() );
    }
  }
}
