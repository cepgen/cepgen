/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2024  Laurent Forthomme
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

#ifndef CepGen_Utils_ProgressBar_h
#define CepGen_Utils_ProgressBar_h

#include <memory>
#include <string>

namespace cepgen::utils {
  class Timer;
  /// A simple progress indicator
  class ProgressBar {
  public:
    /// Progress bar constructor
    /// \param[in] total Total number of steps before completion
    /// \param[in] period Periodicity of ticks refreshing
    explicit ProgressBar(size_t total, size_t period = 10);
    ~ProgressBar();

    /// Reset the progress bar to its initial state
    void reset();
    /// Enable the timer?
    inline void setTimer(bool timer_enabled = true) { timer_enabled_ = timer_enabled; }
    /// Broadcast the current progress to the bar
    /// \param[in] iter Current iteration
    void update(size_t iter) const;

  private:
    std::unique_ptr<Timer> timer_;   ///< Time tracker
    const size_t total_;             ///< Total number of iterations expected
    const size_t period_;            ///< Period at which an iteration is reported
    const size_t bar_length_;        ///< Total number of ticks
    const std::string bar_pattern_;  ///< Characters to use for the bar
    const bool enabled_;             ///< Is the progress bar enabled?

    bool timer_enabled_{false};  ///< Do we also track the time?
    mutable bool ended_{false};  ///< Has the counting stopped?
  };
}  // namespace cepgen::utils

#endif
