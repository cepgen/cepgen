/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2024  Laurent Forthomme
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

#ifndef CepGen_Process_Central2to4PhaseSpaceGenerator_h
#define CepGen_Process_Central2to4PhaseSpaceGenerator_h

#include "CepGen/Process/CentralPhaseSpaceGenerator.h"

namespace cepgen {
  /// A 2-to-4 (or 2-to-2 central) phase space generator
  class Central2to4PhaseSpaceGenerator : public CentralPhaseSpaceGenerator {
  public:
    explicit Central2to4PhaseSpaceGenerator(const ParametersList&);

    size_t ndim() const override { return 4; }
    const pdgids_t& particles() const override { return particles_; }

    void initialise() override;
    double generateKinematics() override;

  protected:
    // mapped variables
    double m_y_c1_{0.};         ///< Rapidity of the first central particle
    double m_y_c2_{0.};         ///< Rapidity of the second central particle
    double m_pt_diff_{0.};      ///< Transverse momentum difference for the two central particle
    double m_phi_pt_diff_{0.};  ///< Azimuthal angle difference for the two central particles

  private:
    // factor 1/4 from jacobian of transformations
    static constexpr double prefactor_ = 0.25 * 0.0625 * M_1_PI * M_1_PI;
    const std::vector<int> int_particles_;
    const pdgids_t particles_;  ///< Type of particles produced in the final state
  };
}  // namespace cepgen

#endif
