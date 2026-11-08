/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#include "CepGen/Utils/Math.h"

namespace cepgen::utils {
  template <typename T>
  bool positive(const T& val) {
    return val > T{} && std::isfinite(val);
  }
  template bool positive<double>(const double&);
  template bool positive<float>(const float&);
  template bool positive<int>(const int&);

  double fastHypot(double x, double y) {
    x *= x, y *= y;
    return std::sqrt(x + y);
  }

  double fastHypot(double x, double y, double z) {
    x *= x, y *= y, z *= z;
    return std::sqrt(x + y + z);
  }

  double fastSqrtSqDiff(double x, double y, Normalise normalise) {
    if (std::fabs(x) == std::fabs(y))
      return 0.;
    if (const auto square = (x + y) * (x - y); normalise == Normalise::yes && square < 0.)
      return -std::sqrt(-square);
    else
      return std::sqrt(square);
  }

  double fastSqrtSqDiff(double x, double y, double z, Normalise normalise) {
    if (std::fabs(x) == std::fabs(y))
      return 0.;
    z *= z;
    if (const auto square = (x + y) * (x - y) - z; normalise == Normalise::yes && square < 0.)
      return -std::sqrt(-square);
    else
      return std::sqrt(square);
  }
}  // namespace cepgen::utils
