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

#ifndef CepGenPythia6_DiffractiveState_h
#define CepGenPythia6_DiffractiveState_h

#include <CepGen/Physics/Hadroniser.h>

namespace cepgen::utils {
  class RandomGenerator;
}  // namespace cepgen::utils

namespace cepgen::pythia6 {
  /// Colour-singlet (anti)quark-(anti)diquark string replacing a dissociated beam, in its own rest frame
  struct DiffractiveState final : hadr::Hadroniser::DiffractiveState {
    /// Build the partonic content of a dissociated system
    /// \param[in] rng Random number generator
    /// \param[in] beam Incoming beam particle
    /// \param[in] system Dissociated system momentum
    explicit DiffractiveState(utils::RandomGenerator& rng,
                              const Particle& beam,
                              const Momentum& system,
                              double primordial_kt);

    /// Run the Lund string fragmentation of the diffractive system
    /// \param[out] fragments Hadrons produced, in the rest frame of the system
    /// \param[in] tolerance Relative tolerance on the four-momentum balance of the final state
    /// \return Boolean stating whether the system was successfully hadronised
    bool hadronise(hadr::Hadroniser::Fragments& fragments, double tolerance) const override;
    double massFromPdgId(long pdg_id) const override;
  };
}  // namespace cepgen::pythia6

#endif
