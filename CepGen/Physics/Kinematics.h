/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#include <iosfwd>
#include <memory>
#include <vector>

#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/IncomingBeams.h"

namespace cepgen {
  /// List of kinematic constraints to apply on the process phase space.
  class Kinematics final : public SteeredObject<Kinematics> {
  public:
    Kinematics();
    explicit Kinematics(const ParametersList&);
    ~Kinematics() = default;

    /// Minimal diffractive mass for dissociative proton treatment
    static const double MX_MIN;

    static ParametersDescription description();

    void setParameters(const ParametersList&) override;
    /// List containing all parameters handled
    ParametersList allParameters(bool extended) const;

    /// Beam/primary particle's kinematics
    IncomingBeams& incomingBeams() { return incoming_beams_; }
    /// Const-qualified beam/primary particle's kinematics
    const IncomingBeams& incomingBeams() const { return incoming_beams_; }

    typedef std::vector<pdgid_t> pdgids_t;
    /// Minimum list of central particles required
    const pdgids_t& minimumFinalState() const { return minimum_final_state_; }

    /// Phase space cuts
    CutsList& cuts() { return cuts_; }
    /// Const-qualified phase space cuts
    const CutsList& cuts() const { return cuts_; }

  private:
    /// Beam/primary particle's kinematics
    IncomingBeams incoming_beams_;
    CutsList cuts_;  ///< Phase space cuts
    pdgids_t minimum_final_state_;
  };
}  // namespace cepgen

#endif
