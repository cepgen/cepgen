/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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

#include <cmath>
#include <sstream>

#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

using namespace cepgen::utils;

void TimeKeeper::clear() {
  monitors_.clear();
  tmr_.reset();
}

TimeKeeper& TimeKeeper::tick(const std::string& func, double time) {
  monitors_[func].emplace_back(time > 0. ? time : tmr_.elapsed());
  return *this;
}

std::string TimeKeeper::summary() const {
  if (monitors_.empty())
    return {};

  // a bit of arithmetics
  struct Monitor {
    std::string name;
    size_t size;
    double total, mean, rms;
    bool operator<(const Monitor& oth) const { return total < oth.total; }
  };

  std::vector<Monitor> mons;
  double total_time = 0.;
  for (const auto& mon : monitors_) {
    const auto& tm = mon.second;
    const double total = tm.empty() ? -1. : std::accumulate(tm.begin(), tm.end(), 0.);
    const double mean = total / static_cast<double>(tm.size());
    const double rms = std::sqrt(std::fabs(
        std::inner_product(tm.begin(), tm.end(), tm.begin(), 0.) / static_cast<double>(tm.size()) - mean * mean));
    mons.emplace_back(Monitor{mon.first, tm.size(), total, mean, rms});
    total_time += total;
  }

  std::sort(mons.rbegin(), mons.rend());  // sort by total clock time (desc.)

  // display the various probes
  static constexpr double s_to_ms = 1.e3;
  std::ostringstream oss;
  oss << format("%2s | %-100s | %12s\t%10s\t%5s", "#", "Caller", "Total (ms)", "Average (ms)", "RMS (ms)");
  for (const auto& mon : mons)
    oss << format("\n%10u | %-100s | %12.6f\t%10e\t%5.3e",
                  mon.size,
                  mon.name.c_str(),
                  mon.total * s_to_ms,
                  mon.mean * s_to_ms,
                  mon.rms * s_to_ms);
  oss << "\nTotal time: " << total_time << ".";
  return oss.str();
}

TimeKeeper::Ticker::Ticker(TimeKeeper* tk, const std::string& name) : tk_(tk), name_(name) {}

TimeKeeper::Ticker::~Ticker() {
  if (tk_)
    tk_->tick(name_, tmr_.elapsed());
}
