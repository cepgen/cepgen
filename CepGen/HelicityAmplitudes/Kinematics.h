/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#ifndef CepGen_HelicityAmplitudes_Kinematics_h
#define CepGen_HelicityAmplitudes_Kinematics_h

#include <array>
#include <iosfwd>

namespace cepgen::amplitudes {
  /// Helper container to handle all kinematics variables computation once
  class Kinematics {
  public:
    explicit Kinematics(double mv2, double s_hat, double t_hat, double u_hat);

    /// Build a kinematics block from its angular and momenta information
    static Kinematics fromSCosTheta(double shat, double cos_theta, double mv2);

    bool operator!=(const Kinematics&) const;
    friend std::ostream& operator<<(std::ostream&, const Kinematics&);

    // base variables
    const double shat{0.}, that{0.}, uhat{0.};

  private:
    void setCosTheta(double cos_theta);  ///< Set the rotation angle (in rad)
    const double mv2_{0.};               /// Squared mass of each of the two produced bosons (in GeV^2)

  public:
    // all derived variables
    const double shat2{0.};
    const double beta2{0.}, beta{0.};
    const double inv_gamma2{0.}, gamma2{0.}, gamma{0.}, inv_gamma{0.};
    double cos_theta{0.}, cos_theta2{0.}, sin_theta2{0.}, sin_theta{0.};
    double invA{0.};
  };

  using HelicityState = std::array<short, 4>;  ///< Two incoming, and two outgoing branches
}  // namespace cepgen::amplitudes

#endif
