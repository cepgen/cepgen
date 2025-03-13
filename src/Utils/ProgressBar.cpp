/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2025  Laurent Forthomme
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

#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Logger.h"
#include "CepGen/Utils/ProgressBar.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Timer.h"

using namespace cepgen::utils;

ProgressBar::ProgressBar(size_t total, size_t period)
    : timer_(new Timer),
      total_(total),
      period_(period),
      bar_length_(std::stoi(env::get("COLUMNS", "60")) - 10),
      bar_pattern_(bar_length_, '='),
      enabled_(env::get("CG_CI").empty() && Logger::get().isTTY()) {}

ProgressBar::~ProgressBar() = default;

void ProgressBar::reset() {
  timer_->reset();
  ended_ = false;
}

void ProgressBar::update(size_t iter) const {
  if (!enabled_ || ended_)
    return;
  if (iter >= total_ - period_) {  // counter is over
    const std::string message = format("[Finished in %g s]", timer_->elapsed());
    fflush(stderr);
    fprintf(stderr,
            "\r%s%.*s%*s\n",
            message.data(),
            0,
            "",
            static_cast<int>(bar_length_ + (timer_enabled_ ? extra_bar_length_ : 0)),
            "");
    fflush(stderr);
    ended_ = true;
    return;
  }
  if (const int percent = iter * 100. / total_; percent % period_ == 0) {  // update the bar displayed
    std::string extra_text;
    if (timer_enabled_) {
      const auto elapsed_time = timer_->elapsed(), expected_time = elapsed_time * total_ / iter;
      extra_text =
          " " + format("%.2fs/%.2fs (remaining: %.2fs)", elapsed_time, expected_time, expected_time - elapsed_time);
      extra_text.resize(extra_bar_length_);
    }
    const auto left_padding = static_cast<int>(percent / 100. * bar_length_),
               right_padding = static_cast<int>(bar_length_) - left_padding;
    fprintf(stderr,
            "\r%3i%% [%.*s%*s]%s",
            percent,
            left_padding,
            bar_pattern_.c_str(),
            right_padding,
            "",
            extra_text.c_str());
    fflush(stderr);
  }
}
