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

#include <string>

namespace cepgen::utils {
  class Timer;
  /// A simple progress indicator
  class ProgressBar {
  public:
    /// Progress bar constructor
    /// \param[in] tot Total number of steps before completion
    /// \param[in] freq Frequency at which the tick is refreshed
    explicit ProgressBar(size_t tot, size_t freq = 10);
    ~ProgressBar();
    /// Broadcast the current progress to the bar
    /// \param[in] iter Current iteration
    void update(size_t iter) const;

  private:
    std::unique_ptr<Timer> tmr_;
    const size_t bar_length_;
    const std::string bar_pattern_;
    const bool enabled_;
    size_t total_, frequency_;
  };
}  // namespace cepgen::utils

#endif
