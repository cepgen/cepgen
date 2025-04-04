/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2014-2025  Laurent Forthomme
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

#ifndef CepGen_Physics_Kinematics_h
#define CepGen_Physics_Kinematics_h

#include "CepGen/Physics/CutsList.h"
#include "CepGen/Physics/IncomingBeams.h"

namespace cepgen {
  /// List of kinematic constraints to apply on the process phase space.
  class Kinematics final : public SteeredObject<Kinematics> {
  public:
    explicit Kinematics(const ParametersList&);

    static ParametersDescription description();

    void setParameters(const ParametersList&) override;
    const ParametersList& parameters() const override;  ///< List containing all parameters handled

    inline IncomingBeams& incomingBeams() { return incoming_beams_; }              ///< Beam/primary particle kinematics
    inline const IncomingBeams& incomingBeams() const { return incoming_beams_; }  ///< Beam/primary particle kinematics

    inline const pdgids_t& minimumFinalState() const { return minimum_final_state_; }  ///< Minimum central particles

    inline CutsList& cuts() { return cuts_; }              ///< Phase space cuts
    inline const CutsList& cuts() const { return cuts_; }  ///< Phase space cuts

  private:
    IncomingBeams incoming_beams_{ParametersList()};  ///< Beam/primary particle's kinematics
    CutsList cuts_{ParametersList()};                 ///< Phase space cuts
    pdgids_t minimum_final_state_;                    ///< Minimum list of particle ids to find in the final state
  };
}  // namespace cepgen

#endif
