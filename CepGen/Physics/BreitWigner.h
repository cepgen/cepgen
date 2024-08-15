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
    explicit BreitWigner(double mean = 0., double gamma = 0., double min_energy = -1., double max_energy = -1.);

    double min() const { return min_energy_; }  ///< Minimal energy to consider
    double max() const { return max_energy_; }  ///< Maximal energy to consider

    double operator()(double x) const;  ///< Shoot a value according to parameterisation

  private:
    const double mean_;        ///< Mean of distribution
    const double gamma_;       ///< Width of distribution
    const double min_energy_;  ///< Minimum energy
    const double max_energy_;  ///< Maximum energy
  };
}  // namespace cepgen

#endif
