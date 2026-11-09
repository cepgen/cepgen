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

#ifndef CepGenHerwig6_Fragment_h
#define CepGenHerwig6_Fragment_h

#include <CepGen/Physics/Momentum.h>

#include <vector>

namespace cepgen::herwig6 {
  struct DiffractiveSystemSeed;

  /// Hadron retrieved from the Herwig 6 event record
  struct Fragment {
    int index{0};        ///< Position in the event record
    int mother{-1};      ///< Position of the first mother in the event record
    long pdg_id{0};      ///< Signed PDG identifier
    bool stable{false};  ///< Is this particle in the final state?
    Momentum momentum;   ///< 4-momentum in the rest frame of the dissociated system
  };
  using Fragments = std::vector<Fragment>;

  /// Run cluster hadronisation of a given diffractive system
  /// \param[in] seed Partonic content of the system
  /// \param[out] fragments Hadrons produced, in the rest frame of the system
  /// \return Boolean stating whether the event was successfully hadronised
  bool hadronise(const DiffractiveSystemSeed& seed, Fragments& fragments);
}  // namespace cepgen::herwig6

#endif
