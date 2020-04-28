#ifndef CepGen_Utils_TimeKeeper_h
#define CepGen_Utils_TimeKeeper_h

#include "CepGen/Utils/Timer.h"

#include <unordered_map>
#include <vector>

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
        TimeKeeper& tick( const std::string& func );
        /// Write a summary of all monitors
        std::string summary() const;

      private:
        std::unordered_map<std::string,std::vector<float> > timers_;
        Timer tmr_;
    };
  }
}

#endif
