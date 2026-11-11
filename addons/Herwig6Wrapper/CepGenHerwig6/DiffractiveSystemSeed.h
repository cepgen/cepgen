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

#ifndef CepGenHerwig6_DiffractiveSystemSeed_h
#define CepGenHerwig6_DiffractiveSystemSeed_h

#include <CepGen/Physics/Hadroniser.h>
#include <CepGen/Physics/Momentum.h>

#include <array>

namespace cepgen::utils {
  class RandomGenerator;
}

namespace cepgen::herwig6 {
  /// Colour-singlet (anti)quark-(anti)diquark system replacing a dissociated beam, in its own rest frame
  class DiffractiveSystemSeed final : public hadr::Hadroniser::DiffractiveState {
  public:
    /// Build the partonic content of a dissociated system
    /// \param[in] beam Incoming beam particle momentum
    /// \param[in] system Dissociated system momentum
    /// \param[in] beam_pdg_id Signed PDG identifier of the incoming beam (only protons and antiprotons are supported)
    /// \param[in] rmass Herwig 6 constituent masses, indexed by Herwig identifier
    /// \param[in] rng Random number generator
    /// \param[in] primordial_kt RMS of the transverse momentum of the spectator relative to the beam axis, in GeV
    static DiffractiveSystemSeed build(const Momentum& beam,
                                       const Momentum& system,
                                       long beam_pdg_id,
                                       const double* rmass,
                                       utils::RandomGenerator& rng,
                                       double primordial_kt);

    /// Run cluster hadronisation of a given diffractive system
    /// \param[out] fragments Hadrons produced, in the rest frame of the system
    /// \param[in] tolerance Relative tolerance on the four-momentum balance of the final state
    /// \return Boolean stating whether the event was successfully hadronised
    bool hadronise(hadr::Hadroniser::Fragments& fragments, double tolerance) const override;

    /// Herwig 6 internal identifiers used for the beam remnants treatment
    enum struct InternalId {
      invalid = 0,
      down = 1,
      up = 2,
      dBar = 7,
      uBar = 8,
      centreOfMass = 14,
      piZero = 21,
      proton = 73,
      uu = 109,
      ud0 = 110,
      uu_bar = 115,
      ud0_bar = 116
    };

  private:
    InternalId struck_id_{InternalId::invalid};     ///< Herwig identifier of the struck (anti)quark
    InternalId spectator_id_{InternalId::invalid};  ///< Herwig identifier of the spectator (anti)diquark
  };
}  // namespace cepgen::herwig6

#endif
