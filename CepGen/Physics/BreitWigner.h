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

#ifndef CepGen_Physics_BreitWigner_h
#define CepGen_Physics_BreitWigner_h

namespace cepgen {
  /// A Breit-Wigner/Cauchy distribution generator
  class BreitWigner {
  public:
    explicit BreitWigner(double mean = 0., double gamma = 0., double emin = -1., double emax = -1.);

    /// Minimal energy to consider
    double min() const { return emin_; }
    /// Maximal energy to consider
    double max() const { return emax_; }
    /// Shoot a value according to parameterisation
    double operator()(double x) const;

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
