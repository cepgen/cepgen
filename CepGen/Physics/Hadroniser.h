/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

    /// Hadron retrieved from the Herwig 6 event record
    struct Fragment {
      int index{0};        ///< Position in the event record
      int mother{-1};      ///< Position of the first mother in the event record
      long pdg_id{0};      ///< Signed PDG identifier
      bool stable{false};  ///< Is this particle in the final state?
      Momentum momentum;   ///< 4-momentum in the rest frame of the dissociated system
    };
    using Fragments = std::vector<Fragment>;

    inline bool fragmentRemnants() const { return fragment_remnants_; }  ///< Fragment the beam remnants?

  protected:
    const bool fragment_remnants_;  ///< Switch on/off the remnants fragmentation where applicable
  };
}  // namespace cepgen::hadr

#endif
