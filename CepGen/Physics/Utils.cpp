/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include "CepGen/Physics/Utils.h"
#include "CepGen/Utils/Limits.h"
#include "CepGen/Utils/Math.h"

namespace cepgen {
  namespace utils {
    static const Limits x_limits{0., 1.};

    double mX2(double xbj, double q2, double mp2) {
      if (!x_limits.contains(xbj))
        return 0.;
      return mp2 + q2 * (1. - xbj) / xbj;
    }

    double q2(double xbj, double mp2, double mx2) {
      if (!x_limits.contains(xbj))
        return 0.;
      return xbj / (1. - xbj) * (mx2 - mp2);
    }

    double xBj(double q2, double mp2, double mx2) {
      if (!positive(q2))
        return 0.;
      return q2 / (q2 + mx2 - mp2);
    }

    double energyFromW(double w, double mp2, double m2) {
      if (!positive(w))
        return 0.;
      return 0.5 * (w * w - mp2 + m2) / w;
    }

    namespace kt {
      double mX2(double x, double kt2, double q2, double mi2) {
        if (!positive(x))
          return 0.;
        return mi2 + (q2 * (1. - x) - kt2 - x * x * mi2) / x;
      }

      double q2(double x, double kt2, double mi2, double mx2) {
        if (!x_limits.contains(x))
          return 0.;
        if (mx2 < 0.)  // mx2 = mi2
          return (kt2 + x * x * mi2) / (1. - x);
        return (kt2 + x * (mx2 - mi2) + x * x * mi2) / (1. - x);
      }
    }  // namespace kt
  }    // namespace utils
}  // namespace cepgen
