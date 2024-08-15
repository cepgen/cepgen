/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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

#include "CepGen/Utils/Histogram.h"

#include "CepGen/Utils/Message.h"

namespace cepgen::utils {
  void Histogram::normalise(double integral_value) { scale(integral_value / integral()); }

  std::set<double> Histogram::extractBins(BinMode mode,
                                          size_t num_bins,
                                          const std::function<Limits(size_t)>& bins_extractor) const {
    std::set<double> bins;
    for (size_t i = 0; i < num_bins; ++i) {
      const auto range = bins_extractor(i);
      if (mode == BinMode::low || mode == BinMode::both)
        bins.insert(range.min());
      if (mode == BinMode::high || mode == BinMode::both)
        bins.insert(range.max());
    }
    const auto exp_bins = mode == BinMode::both ? num_bins + 1 : num_bins;
    if (bins.size() != exp_bins)
      CG_WARNING("Histogram:extractBins")
          << "Invalid number of values to bin ranges: got " << bins.size() << ", expecting " << exp_bins << ".";
    return bins;
  }
}  // namespace cepgen::utils
