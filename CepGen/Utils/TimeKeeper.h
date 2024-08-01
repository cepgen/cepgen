/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2024  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Utils_TimeKeeper_h
#define CepGen_Utils_TimeKeeper_h

#include <string>
#include <unordered_map>
#include <vector>

#include "CepGen/Utils/Timer.h"

#define CG_CONCAT(a, b) a##b
#define CG_TICKER_NAME(a, b) CG_CONCAT(a, b)
#define CG_TICKER(tmr) utils::TimeKeeper::Ticker CG_TICKER_NAME(ticker, __COUNTER__)(tmr, __PRETTY_FUNCTION__)

namespace cepgen::utils {
  /// Collection of clocks to benchmark execution blocks
  class TimeKeeper {
  public:
    explicit TimeKeeper() = default;

    void clear();                                            ///< Reset all counters and the timer
    inline bool empty() const { return monitors_.empty(); }  ///< Check if at least one monitor recorded something
    std::string summary() const;                             ///< Write a summary of all monitors

    /// Count the time for one monitor
    /// \param[in] func monitor to increment
    /// \param[in] time increment, in second (< 0 to count since last timer reset)
    TimeKeeper& tick(const std::string& func, double time = -1.);

    const Timer& timer() const;  ///< Local timer object

    /// Scoped timekeeping utility
    class Ticker {
    public:
      explicit Ticker(TimeKeeper*, const std::string&);  ///< Build a named and scoped time ticker
      ~Ticker();  ///< Ticker destructor to store the timing information to the parent timekeeper

    private:
      TimeKeeper* tk_{nullptr};
      std::string name_;
      Timer tmr_;
    };

  private:
    std::unordered_map<std::string, std::vector<float> > monitors_;
    Timer tmr_;
  };
}  // namespace cepgen::utils

#endif
