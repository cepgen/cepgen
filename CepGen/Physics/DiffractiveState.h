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

#ifndef CepGen_Physics_DiffractiveState_h
#define CepGen_Physics_DiffractiveState_h

#include <vector>

#include "CepGen/Physics/Momentum.h"

namespace cepgen::utils {
  class RandomGenerator;
}  // namespace cepgen::utils

/// Location for all hadronisers to be run downstream to the events generation
namespace cepgen::hadr {
  /// Hadron retrieved from the hadroniser event record
  struct Fragment {
    int index{0};        ///< Position in the event record
    int mother{-1};      ///< Position of the first mother in the event record
    long pdg_id{0};      ///< Signed PDG identifier
    bool stable{false};  ///< Is this particle in the final state?
    Momentum momentum;   ///< 4-momentum in the rest frame of the dissociated system
  };
  using Fragments = std::vector<Fragment>;

  /// Colour-singlet (anti)quark-(anti)diquark system replacing a dissociated beam, in its own rest frame
  class DiffractiveState {
  public:
    struct ContentProbability {
      long valquark_pdg_id{0};  ///< PDG identifier of the struck quark
      long diquark_pdg_id{0};   ///< PDG identifier of the spectator diquark
      double weight{0.};        ///< Relative probability
    };
    using PartonsContent = std::vector<ContentProbability>;

    explicit DiffractiveState(utils::RandomGenerator& rng, long beam_pdg_id, const PartonsContent&);
    virtual ~DiffractiveState() = default;

    /// Run the diffractive system fragmentation
    /// \param[out] fragments Hadrons produced, in the rest frame of the system
    /// \param[in] tolerance Relative tolerance on the four-momentum balance of the final state
    /// \return Boolean stating whether the diffractive state was successfully hadronised
    virtual bool hadronise(Fragments& fragments, double tolerance) const = 0;

    bool valid() const { return valid_; }  ///< Is the fragmentation of the diffractive state possible?
  protected:
    /// Compute the parton masses and the system kinematics
    void initialise(const Momentum& system_momentum, const Momentum& beam_momentum, double primordial_kt);
    /// Check the momentum balance between the system and its stable fragments
    bool momentumConserved(const Fragments& fragments, double tolerance) const;
    /// Mass of a particle according to the hadroniser algorithm
    virtual double massFromPdgId(long pdg_id) const;

    utils::RandomGenerator* random_generator_{nullptr};  //NOT owning

    bool valid_{true};            ///< Is the fragmentation of the diffractive state possible?
    long beam_pdg_id_{0};         ///< Signed PDG identifier of the dissociated beam
    long valquark_pdg_id_{0};     ///< Signed PDG identifier of the struck (anti)quark
    long diquark_pdg_id_{0};      ///< Signed PDG identifier of the spectator (anti)diquark
    double system_mass_{0.};      ///< Invariant mass of the system
    double valquark_mass_{0.};    ///< Mass of the struck (anti)quark
    double diquark_mass_{0.};     ///< Mass of the spectator (anti)diquark
    size_t mass_model_index_{0};  ///< Index of the mass model eventually used
    Momentum diquark_momentum_;   ///< Spectator 3-momentum (struck parton has the opposite)
  };
}  // namespace cepgen::hadr

#endif
