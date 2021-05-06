#ifndef CepGen_Utils_TimeKeeper_h
#define CepGen_Utils_TimeKeeper_h

#include "CepGen/Utils/Timer.h"

#include <string>
#include <unordered_map>
#include <vector>

#define CG_CONCAT(a, b) a##b
#define CG_TICKER_NAME(a, b) CG_CONCAT(a, b)
#define CG_TICKER(tmr) utils::TimeKeeper::Ticker CG_TICKER_NAME(ticker, __COUNTER__)(tmr, __PRETTY_FUNCTION__)

namespace cepgen {
  namespace utils {
    /// A collection of clocks to benchmark execution blocks
    class TimeKeeper {
    public:
      /// Object constructor
      explicit TimeKeeper();

      /// Reset all counters and the timer
      void clear();
      /// Check if at least one monitor recorded something
      bool empty() const { return monitors_.empty(); }
      /// Count the time for one monitor
      /// \param[in] func Which monitor to increment
      /// \param[in] time Increment, in second (< 0 to count since last timer reset)
      TimeKeeper& tick(const std::string& func, double time = -1.);
      /// Write a summary of all monitors
      std::string summary() const;

      /// Local timer object
      const Timer& timer() const;
      /// A scoped timekeeping utility
      class Ticker {
      public:
        /// Build a named and scoped time ticker
        explicit Ticker(TimeKeeper*, const std::string&);
        /// Ticker destructor to store the timing information to the parent timekeeper
        ~Ticker();

      private:
        TimeKeeper* tk_;
        std::string name_;
        Timer tmr_;
      };

    private:
      std::unordered_map<std::string, std::vector<float> > monitors_;
      Timer tmr_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
