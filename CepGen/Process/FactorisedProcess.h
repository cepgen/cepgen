/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#ifndef CepGen_Process_FactorisedProcess_h
#define CepGen_Process_FactorisedProcess_h

#include "CepGen/Process/PhaseSpaceGenerator.h"
#include "CepGen/Process/Process.h"

namespace cepgen {
  class PartonFlux;
  namespace proc {
    /**
     * A generic factorised process.
     * \note 0 to 2 dimensions may be used for the scattered diffractive system(s)' invariant mass definition.
     * \brief Class template to define any factorised process
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2023
     */
    class FactorisedProcess : public Process {
    public:
      /// Class constructor
      /// \param[in] params Parameters list
      /// \param[in] output Produced final state particles
      explicit FactorisedProcess(const ParametersList& params, const pdgids_t& output);
      FactorisedProcess(const FactorisedProcess&);

      double computeWeight() override final;
      void fillKinematics(bool) override final;

      static ParametersDescription description();

    protected:
      void addEventContent() override final;
      void prepareKinematics() override final;

      /// Set the kinematics of the central system before any point computation
      virtual void setExtraContent() {}
      /// Prepare the central part of the Jacobian (only done once, as soon as the kinematics is set)
      virtual void preparePhaseSpace() = 0;
      /// Factorised matrix element (event weight)
      virtual double computeFactorisedMatrixElement() = 0;
      /// Set the kinematics of the outgoing central system
      virtual void fillCentralParticlesKinematics() = 0;

      /// Type of particles produced in the final state
      pdgids_t produced_parts_;
      /// Kinematic variables generator for the phase space coverage
      const std::unique_ptr<PhaseSpaceGenerator> psgen_;
    };
  }  // namespace proc
}  // namespace cepgen

#endif
