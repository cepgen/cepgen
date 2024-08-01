/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

#ifndef CepGen_Utils_Math_h
#define CepGen_Utils_Math_h

namespace cepgen::utils {
  /// Check if a number is positive and finite
  template <typename T>
  bool positive(const T& val);

  double fastHypot(double, double);
  double fastHypot(double, double, double);
  /// Compute the square root of the squared difference (sqrt(a^2-b^2))
  double fastSqrtSqDiff(double, double);
}  // namespace cepgen::utils

#endif
