/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2026  Laurent Forthomme
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

#ifndef CepGen_Physics_Hadroniser_h
#define CepGen_Physics_Hadroniser_h

#include <array>

#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Physics/Momentum.h"

/// Location for all hadronisers to be run downstream to the events generation
namespace cepgen::hadr {
  /// Class template to define any hadroniser as a general object with defined methods
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date January 2014
  class Hadroniser : public EventModifier {
  public:
    explicit Hadroniser(const ParametersList&);  ///< Default constructor for an undefined hadroniser

    static ParametersDescription description();

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
      explicit DiffractiveState(const Momentum& system_momentum);

      /// Run the diffractive system fragmentation
      /// \param[out] fragments Hadrons produced, in the rest frame of the system
      /// \param[in] tolerance Relative tolerance on the four-momentum balance of the final state
      /// \return Boolean stating whether the diffractive state was successfully hadronised
      virtual bool hadronise(hadr::Hadroniser::Fragments& fragments, double tolerance) const = 0;

      bool valid() const { return valid_; }  ///< Is the fragmentation of the diffractive state possible?

    protected:
      bool valid_{true};                                    ///< Is the fragmentation of the diffractive state possible?
      long valquark_pdg_id_{0};                             ///< Signed PDG identifier of the struck (anti)quark
      long diquark_pdg_id_{0};                              ///< Signed PDG identifier of the spectator (anti)diquark
      double system_mass_{0.};                              ///< Invariant mass of the system
      double valquark_mass_{0.};                            ///< Mass of the struck (anti)quark
      double diquark_mass_{0.};                             ///< Mass of the spectator (anti)diquark
      std::array<double, 3> diquark_momentum_{0., 0., 0.};  ///< Spectator 3-momentum (struck parton has the opposite)
                                                            //TODO: use Momentum here
    };

    inline bool fragmentRemnants() const { return fragment_remnants_; }  ///< Fragment the beam remnants?

  protected:
    const bool fragment_remnants_;  ///< Switch on/off the remnants fragmentation where applicable
  };
}  // namespace cepgen::hadr

#endif
