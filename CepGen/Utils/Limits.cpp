#include "CepGen/Utils/Limits.h"
#include "CepGen/Utils/String.h"

#include "CepGen/Core/Exception.h"

namespace cepgen {
  Limits::Limits(double min, double max) : std::pair<double, double>(min, max) {}

  Limits::Limits(const Limits& rhs) : std::pair<double, double>(rhs.first, rhs.second) {}

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

  bool Limits::contains(double val) const {
    if (hasMin() && val < min())
      return false;
    if (hasMax() && val > max())
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

  void Limits::save(bool& on, double& lmin, double& lmax) const {
    on = false;
    lmin = lmax = 0.;
    if (!valid())
      return;
    on = true;
    if (hasMin())
      lmin = min();
    if (hasMax())
      lmax = max();
    if (lmin == lmax)
      on = false;
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
