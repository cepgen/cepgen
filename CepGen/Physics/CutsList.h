/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2024  Laurent Forthomme
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

#ifndef CepGen_Physics_CutsList_h
#define CepGen_Physics_CutsList_h

#include "CepGen/Physics/Cuts.h"

namespace cepgen {
  /// Collection of cuts to be applied on all particle with a given PDG id
  using PerIdCuts = std::unordered_map<pdgid_t, cuts::Central>;
  /// A collection of cuts to apply on the physical phase space
  struct CutsList final : SteeredObject<CutsList> {
    explicit CutsList(const ParametersList&);

    static ParametersDescription description();

    void setParameters(const ParametersList&) override;

    friend std::ostream& operator<<(std::ostream&, const CutsList&);

    cuts::Initial initial;        ///< Cuts on the initial particles kinematics
    cuts::Central central;        ///< Cuts on the central system produced
    PerIdCuts central_particles;  ///< Cuts on the central individual particles
    cuts::Remnants remnants;      ///< Cuts on the beam remnants system
  };
}  // namespace cepgen

#endif
