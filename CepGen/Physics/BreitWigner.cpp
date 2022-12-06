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

#include <math.h>

#include "CepGen/Physics/BreitWigner.h"

namespace cepgen {
  BreitWigner::BreitWigner(double mean, double gamma, double emin, double emax)
      : mean_(mean), gamma_(gamma), emin_(emin), emax_(emax) {}

  double BreitWigner::operator()(double x) const {
    const double val = mean_ + 0.5 * gamma_ * tan((2. * x - 1.) * M_PI_2);
    if (emin_ >= 0. && val < emin_)
      return -1.;
    if (emax_ >= 0. && val > emax_)
      return -1.;
    return val;
  }

}  // namespace cepgen
