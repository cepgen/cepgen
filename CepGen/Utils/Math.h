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

#ifndef CepGen_Utils_Math_h
#define CepGen_Utils_Math_h

namespace cepgen::utils {
  /// Check if a number is positive and finite
  template <typename T>
  bool positive(const T& val);
  /// Type-safe C++ sign function
  template <typename T>
  inline short sign(const T& val) {
    return (T(0) < val) - (val < T(0));
  }
  /// Compute the square root of the squared sum (sqrt(a^2+b^2))
  double fastHypot(double, double);
  /// Compute the square root of the squared sum (sqrt(a^2+b^2+c^2))
  double fastHypot(double, double, double);

  enum struct Normalise { no = 0, yes = 1 };
  /// Compute the square root of the squared difference (sqrt(a^2-b^2))
  double fastSqrtSqDiff(double, double, Normalise normalise = Normalise::no);
  /// Compute the square root of the squared difference (sqrt(a^2-b^2-c^2))
  double fastSqrtSqDiff(double, double, double, Normalise normalise = Normalise::no);
}  // namespace cepgen::utils

#endif
