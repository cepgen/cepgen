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
    if (v < 0. || v > 1.)
      throw CG_ERROR("Limits:shoot") << "x must be comprised between 0 and 1; x value = " << v << ".";

    if (!hasMin() && hasMax())  // if no lower limit, assume 0
      return second * v;

    if (!valid())
      return INVALID;

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

  namespace utils {
    Limits log(const Limits& lim) {
      return Limits{lim.hasMin() ? std::log(lim.min()) : Limits::INVALID,
                    lim.hasMax() ? std::log(lim.max()) : Limits::INVALID};
    }

    Limits log10(const Limits& lim) {
      return Limits{lim.hasMin() ? std::log10(lim.min()) : Limits::INVALID,
                    lim.hasMax() ? std::log10(lim.max()) : Limits::INVALID};
    }

    Limits pow(const Limits& lim, double exp) {
      return Limits{lim.hasMin() ? std::pow(lim.min(), exp) : Limits::INVALID,
                    lim.hasMax() ? std::pow(lim.max(), exp) : Limits::INVALID};
    }

    Limits sqrt(const Limits& lim) {
      return Limits{lim.hasMin() ? std::sqrt(lim.min()) : Limits::INVALID,
                    lim.hasMax() ? std::sqrt(lim.max()) : Limits::INVALID};
    }
  }  // namespace utils
}  // namespace cepgen
