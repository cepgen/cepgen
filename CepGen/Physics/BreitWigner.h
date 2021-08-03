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

#ifndef CepGen_Physics_BreitWigner_h
#define CepGen_Physics_BreitWigner_h

#include <math.h>

namespace cepgen {
  /// A Breit-Wigner/Cauchy distribution generator
  class BreitWigner {
  public:
    BreitWigner(double mean = 0., double gamma = 0., double emin = -1., double emax = -1.)
        : mean_(mean), gamma_(gamma), emin_(emin), emax_(emax) {}
    /// Copy constructor
    BreitWigner(const BreitWigner& oth) : mean_(oth.mean_), gamma_(oth.gamma_), emin_(oth.emin_), emax_(oth.emax_) {}

    /// Minimal energy to consider
    double min() const { return emin_; }
    /// Maximal energy to consider
    double max() const { return emax_; }
    /// Shoot a value according to parameterisation
    inline double operator()(double x) const {
      const double val = mean_ + 0.5 * gamma_ * tan((2. * x - 1.) * M_PI_2);
      if (emin_ >= 0. && val < emin_)
        return -1.;
      if (emax_ >= 0. && val > emax_)
        return -1.;
      return val;
    }

  private:
    /// Mean of distribution
    double mean_;
    /// Width of distribution
    double gamma_;
    /// Minimal value
    double emin_;
    /// Maximal value
    double emax_;
  };
}  // namespace cepgen

#endif
