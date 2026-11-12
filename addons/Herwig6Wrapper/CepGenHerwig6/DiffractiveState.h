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

#ifndef CepGenHerwig6_DiffractiveState_h
#define CepGenHerwig6_DiffractiveState_h

#include <CepGen/Physics/DiffractiveState.h>

namespace cepgen {
  class Particle;
}  // namespace cepgen
namespace cepgen::utils {
  class RandomGenerator;
}  // namespace cepgen::utils

namespace cepgen::herwig6 {
  /// Colour-singlet (anti)quark-(anti)diquark system replacing a dissociated beam, in its own rest frame
  class DiffractiveState final : public hadr::DiffractiveState {
  public:
    /// Build the partonic content of a dissociated system
    /// \param[in] rng Random number generator
    /// \param[in] beam Incoming beam particle
    /// \param[in] system Dissociated system momentum
    /// \param[in] primordial_kt RMS of the transverse momentum of the spectator relative to the beam axis, in GeV
    explicit DiffractiveState(utils::RandomGenerator& rng,
                              const Particle& beam,
                              const Momentum& system,
                              double primordial_kt);

    /// Run cluster hadronisation of a given diffractive system
    /// \param[out] fragments Hadrons produced, in the rest frame of the system
    /// \param[in] tolerance Relative tolerance on the four-momentum balance of the final state
    /// \return Boolean stating whether the event was successfully hadronised
    bool hadronise(hadr::Fragments& fragments, double tolerance) const override;

    /// Herwig 6 internal identifiers used for the beam remnants treatment
    enum struct InternalId {
      invalid = 0,
      down = 1,
      up = 2,
      down_bar = 7,
      up_bar = 8,
      centreOfMass = 14,
      piZero = 21,
      proton = 73,
      uu = 109,
      ud0 = 110,
      uu_bar = 115,
      ud0_bar = 116
    };

  private:
    double massFromPdgId(long pdg_id) const override;
    InternalId valquark_id_{InternalId::invalid};  ///< Herwig identifier of the struck (anti)quark
    InternalId diquark_id_{InternalId::invalid};   ///< Herwig identifier of the spectator (anti)diquark
  };
}  // namespace cepgen::herwig6

#endif
