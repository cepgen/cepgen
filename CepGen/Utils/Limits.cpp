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

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Limits.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  Limits::Limits(double min, double max) : std::pair<double, double>(min, max) {}

  Limits::Limits(const Limits& rhs) : std::pair<double, double>(rhs.first, rhs.second) {}

  bool Limits::operator<(const Limits& oth) const {
    if (first < oth.first)
      return true;
    return second < oth.second;
  }

  Limits Limits::operator-() const { return Limits(-second, -first); }

  Limits& Limits::operator+=(double c) {
    first += c;
    second += c;
    return *this;
  }

  Limits& Limits::operator-=(double c) {
    first -= c;
    second -= c;
    return *this;
  }

  Limits& Limits::operator*=(double c) {
    first *= c;
    second *= c;
    if (c < 0.)
      std::swap(first, second);
    return *this;
  }

  void Limits::in(double low, double up) {
    first = low;
    second = up;
  }

  double Limits::range() const {
    if (!hasMin() && hasMax())  // if no lower limit, assume 0
      return second;
    if (!hasMin() || !hasMax())
      return 0.;
    return second - first;
  }

  bool Limits::hasMin() const { return first != INVALID; }

  bool Limits::hasMax() const { return second != INVALID; }

  bool Limits::contains(double val, bool exclude_boundaries) const {
    if (hasMin() && (val < min() || (exclude_boundaries && val == min())))
      return false;
    if (hasMax() && (val > max() || (exclude_boundaries && val == max())))
      return false;
    return true;
  }

  bool Limits::valid() const {
    if (min() == max())
      return false;
    return hasMin() || hasMax();
  }

  Limits& Limits::validate() {
    if (second < first)
      second = INVALID;
    return *this;
  }

  double Limits::x(double v) const {
    static const Limits x_limits{0., 1.};
    if (!x_limits.contains(v))
      throw CG_ERROR("Limits:shoot") << "x = " << v << " must be inside " << x_limits << ".";
    if (!hasMin() && hasMax()) {  // if no lower limit, assume 0
      CG_WARNING("Limits:shoot") << "Requested to give a value inside interval while no lower limit is set. "
                                 << "Assuming this latter is equal to 0.";
      return second * v;
    }
    if (!valid()) {
      CG_WARNING("Limits:shoot") << "Requested to give a value inside interval although this latter is invalid.";
      return INVALID;
    }
    return first + (second - first) * v;
  }

  std::vector<double> Limits::generate(size_t num_bins, bool log_scale) const {
    std::vector<double> out;
    const auto min_val = (!log_scale) ? min() : std::log10(min());
    const auto rng = ((!log_scale) ? max() - min() : std::log10(max()) - std::log10(min())) / (num_bins - 1);
    for (size_t i = 0; i < num_bins; ++i)
      out.emplace_back((!log_scale) ? min_val + i * rng : std::pow(10, min_val + i * rng));
    return out;
  }

  std::vector<Limits> Limits::split(size_t num_bins, bool log_scale) const {
    const auto& gen = generate(num_bins, log_scale);
    std::vector<Limits> out;
    for (size_t i = 0; i < gen.size() - 1; ++i)
      out.emplace_back(gen.at(i), gen.at(i + 1));
    return out;
  }

  Limits Limits::truncate(const Limits& ext) const {
    auto out = *this;
    if (ext.hasMin() && (!out.hasMin() || out.min() < ext.min()))
      out.min() = ext.min();
    if (ext.hasMax() && (!out.hasMax() || out.max() > ext.max()))
      out.max() = ext.max();
    return out;
  }

  double Limits::trim(double val) const {
    if (hasMin() && val < min())
      return min();
    if (hasMax() && val > max())
      return max();
    return val;
  }

  Limits& Limits::apply(double (*op)(double)) {
    (*this) = compute(op);
    return *this;
  }

  Limits Limits::compute(double (*op)(double)) const {
    return compute([&op](double ext) { return op(ext); });
  }

  std::ostream& operator<<(std::ostream& os, const Limits& lim) {
    if (!lim.hasMin() && !lim.hasMax())
      return os << "no cuts";
    if (!lim.hasMin())
      return os << utils::format("below %g", lim.max());
    if (!lim.hasMax())
      return os << utils::format("above %g", lim.min());
    return os << utils::format("%g to %g", lim.min(), lim.max());
  }

  Limits operator+(Limits lim, double c) {
    lim += c;
    return lim;
  }

  Limits operator-(Limits lim, double c) {
    lim -= c;
    return lim;
  }

  Limits operator*(Limits lim, double c) {
    lim *= c;
    return lim;
  }
}  // namespace cepgen
