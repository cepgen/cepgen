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

#include <optional>

namespace cepgen::utils {
  class RandomGenerator;
}

namespace cepgen::pythia6 {
  /// Colour-singlet (anti)quark-(anti)diquark string replacing a dissociated beam, in its own rest frame
  class DiffractiveState final : public hadr::Hadroniser::DiffractiveState {
  public:
    /// Build the partonic content of a dissociated system
    /// \param[in] system Dissociated system momentum
    /// \param[in] beam_pdg_id Signed PDG identifier of the incoming beam (a negative value yields antiquark content)
    /// \param[in] rng Random number generator
    explicit DiffractiveState(const Momentum& system, long beam_pdg_id, utils::RandomGenerator& rng);

    /// Run the Lund string fragmentation of the diffractive system
    /// \param[out] fragments Hadrons produced, in the rest frame of the system
    /// \param[in] tolerance Relative tolerance on the four-momentum balance of the final state
    /// \return Boolean stating whether the system was successfully hadronised
    bool hadronise(hadr::Hadroniser::Fragments& fragments, double tolerance) const override;

  private:
    long struck_pdg_id_{0};     ///< Signed PDG identifier of the struck (anti)quark
    long spectator_pdg_id_{0};  ///< Signed PDG identifier of the spectator (anti)diquark
  };
}  // namespace cepgen::pythia6

#endif
