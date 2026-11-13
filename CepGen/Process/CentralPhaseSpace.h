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

#ifndef CepGen_Process_CentralPhaseSpace_h
#define CepGen_Process_CentralPhaseSpace_h

#include <vector>

#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Process/IntegrationVariable.h"

namespace cepgen {
  class Kinematics;
}  // namespace cepgen

namespace cepgen {
  /// Central system phase space mapper: builds the momenta of the N central particles produced in the
  /// two-parton collision, given the two-parton system transverse momentum
  /// \note The longitudinal kinematics of the incoming partons (momentum fractions) is not handled here: it is
  ///  deduced by the consumer process from the central system, once built.
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Nov 2026
  class CentralPhaseSpace : public NamedModule<CentralPhaseSpace> {
  public:
    explicit CentralPhaseSpace(const ParametersList&);

    static ParametersDescription description();

    /// Set the central particles content, and retrieve the kinematics cuts to apply on them
    void initialise(const Kinematics&, const spdgids_t& central_particles);

    virtual size_t multiplicity() const = 0;  ///< Number of central particles handled by this mapper
    /// List of integration variables required by this mapper (to be registered by the consumer process)
    virtual std::vector<IntegrationVariable> integrationVariables() = 0;
    /// Build the central particles' momenta for the current phase space point
    /// \param[in] qt_sum Two-parton system transverse momentum
    /// \param[out] central_momenta Central particles' momenta (only filled for a valid point)
    /// \return Jacobian of the central phase space for this point (0 if the point is to be rejected)
    virtual double generate(const Momentum& qt_sum, std::vector<Momentum>& central_momenta) const = 0;

    inline const spdgids_t& particles() const { return particles_; }  ///< Central particles content

  protected:
    const Kinematics& kinematics() const;  ///< Consumer process' kinematics (incl. its cuts)
    /// Per-particle cuts, as steered for this particle species (unconstrained if none was set)
    inline const cuts::Central& particleCuts(size_t i) const { return particle_cuts_.at(i); }
    inline double mass(size_t i) const { return masses_.at(i); }  ///< Central particle mass

  private:
    const Kinematics* kin_{nullptr};  ///< Consumer process' kinematics (not owning)
    spdgids_t particles_;
    std::vector<double> masses_;
    std::vector<cuts::Central> particle_cuts_;
  };
}  // namespace cepgen

#endif
