/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include "CepGen/Utils/ProgressBar.h"

namespace cepgen {
  namespace utils {
    ProgressBar::ProgressBar(size_t tot, size_t freq)
        : bar_length_(std::stoi(env::get("COLUMNS", "60")) - 10),
          bar_pattern_(bar_length_, '='),
          enabled_(env::get("CG_CI").empty()),
          total_(tot),
          frequency_(freq) {}

    void ProgressBar::update(size_t iter) const {
      if (!enabled_)
        return;
      const size_t percent = iter * 100. / total_;
      if (percent % frequency_ == 0 || iter == total_) {
        int lpad = int(percent / 100. * bar_length_);
        int rpad = bar_length_ - lpad;
        fprintf(stderr, "\r%3zu%% [%.*s%*s]", percent, lpad, bar_pattern_.c_str(), rpad, "");
        fflush(stderr);
        if (iter == total_)
          fprintf(stderr, "\n");
      }
    }
  }  // namespace utils
}  // namespace cepgen
