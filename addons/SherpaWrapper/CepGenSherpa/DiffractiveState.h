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

#ifndef CepGenSherpa_DiffractiveState_h
#define CepGenSherpa_DiffractiveState_h

#include <CepGen/Physics/Hadroniser.h>
#include <CepGen/Physics/Momentum.h>

namespace cepgen {
  class Particle;
}  // namespace cepgen
namespace cepgen::utils {
  class RandomGenerator;
}  // namespace cepgen::utils

namespace cepgen::sherpa {
  /// Colour-singlet (anti)quark-(anti)diquark system replacing a dissociated beam, in its own rest frame
  class DiffractiveState final : public hadr::Hadroniser::DiffractiveState {
  public:
    /// Build the partonic content of a dissociated system
    /// \param[in] rng Random number generator
    /// \param[in] beam Incoming beam particle
    /// \param[in] system Dissociated system momentum
    /// \note Partons are put on their constituent mass shells whenever the system is heavy enough, on their
    ///       perturbative mass shells otherwise
    /// \param[in] primordial_kt RMS of the transverse momentum of the spectator relative to the beam axis, in GeV
    explicit DiffractiveState(utils::RandomGenerator& rng,
                              const Particle& beam,
                              const Momentum& system,
                              double primordial_kt);

    /// Run the Sherpa fragmentation (and hadron decays, if enabled) of a given diffractive system
    /// \param[out] fragments Hadrons produced, in the rest frame of the system
    /// \param[in] tolerance Relative tolerance on the four-momentum balance of the final state
    /// \return Boolean stating whether the event was successfully hadronised
    bool hadronise(hadr::Hadroniser::Fragments& fragments, double tolerance) const override;

  private:
    double massFromPdgId(long pdg_id) const override;
    bool constituent_masses_{false};  ///< Are the partons on their constituent mass shell?
  };
}  // namespace cepgen::sherpa

#endif
