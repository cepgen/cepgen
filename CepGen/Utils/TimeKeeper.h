#ifndef CepGen_Utils_TimeKeeper_h
#define CepGen_Utils_TimeKeeper_h

#include "CepGen/Utils/Timer.h"

#include <unordered_map>
#include <vector>

#define CG_TICKER_NAME ticker ## __PRETTY_FUNCTION__
#define CG_TICKER( tmr ) utils::TimeKeeper::Ticker CG_TICKER_NAME( tmr, __PRETTY_FUNCTION__ )

namespace cepgen
{
  namespace utils
  {
    class TimeKeeper
    {
      public:
        explicit TimeKeeper();

        /// Reset all counters and the timer
        void clear();
        /// Count the time for one monitor
        /// \param[in] func Which monitor to increment
        /// \param[in] time Increment, in second (< 0 to count since last timer reset)
        TimeKeeper& tick( const std::string& func, double time = -1. );
        /// Write a summary of all monitors
        std::string summary() const;

        /// Local timer object
        const Timer& timer() const;

        class Ticker
        {
          public:
            explicit Ticker( TimeKeeper&, const std::string& );
            ~Ticker();

          private:
            TimeKeeper& tk_;
            std::string name_;
            Timer tmr_;
        };

      private:
        std::unordered_map<std::string,std::vector<float> > monitors_;
        Timer tmr_;
    };
  }
}

#endif
