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

#include <cmath>

#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Logger.h"
#include "CepGen/Utils/ProgressBar.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Timer.h"

namespace cepgen::utils {
  ProgressBar::ProgressBar(size_t tot, size_t freq)
      : tmr_(new Timer),
        bar_length_(std::stoi(env::get("COLUMNS", "60")) - 10),
        bar_pattern_(bar_length_, '='),
        enabled_(env::get("CG_CI").empty() && Logger::get().isTTY()),
        total_(tot),
        frequency_(freq) {}

  ProgressBar::~ProgressBar() {
    const std::string message = format("[Finished in %g s]", tmr_->elapsed());
    fprintf(stderr, "\r%s%.*s%*s\n", message.data(), 0, "", static_cast<int>(bar_length_), "");
    fflush(stderr);
  }

  void ProgressBar::update(size_t iter) const {
    if (!enabled_)
      return;
    const size_t percent = iter * 100. / total_;
    if (percent % frequency_ == 0 || iter == total_) {
      const auto left_padding = static_cast<int>(percent / 100. * bar_length_),
                 right_padding = static_cast<int>(bar_length_) - left_padding;
      fprintf(stderr, "\r%3zu%% [%.*s%*s]", percent, left_padding, bar_pattern_.c_str(), right_padding, "");
      fflush(stderr);
    }
  }
}  // namespace cepgen::utils
