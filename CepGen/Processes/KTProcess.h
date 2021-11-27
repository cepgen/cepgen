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

#ifndef CepGen_Processes_KTProcess_h
#define CepGen_Processes_KTProcess_h

#include "CepGen/Processes/Process.h"

namespace cepgen {
  namespace proc {
    /**
     * A generic \f$k_{\rm T}\f$-factorisation process.
     * \note
     * - First 4 dimensions of the phase space are required for the
     *    incoming partons' virtualities (radial and azimuthal coordinates).
     * - Last 0-2 dimensions may be used for the scattered diffractive
     *    system(s)' invariant mass definition.
     * \brief Class template to define any \f$k_{\rm T}\f$-factorisation process
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Apr 2016
     */
    class KTProcess : public Process {
    public:
      /// Class constructor
      /// \param[in] params Parameters list
      /// \param[in] partons First and second incoming parton
      /// \param[in] output Produced final state particles
      KTProcess(const ParametersList& params,
                const std::array<pdgid_t, 2>& partons,
                const std::vector<pdgid_t>& output);
      /// Copy constructor
      KTProcess(const KTProcess&);

      /// Populate the event content with the generated process' topology
      void addEventContent() override;
      /// Retrieve the event weight in the phase space
      double computeWeight() override;
      /// Populate the event content with the generated process' kinematics
      void fillKinematics(bool) override;

    protected:
      /// Set the kinematics associated to the phase space definition
      void prepareKinematics() override;
      /// Set the kinematics of the central system before any point computation
      virtual void setExtraContent() {}
      /// Prepare the central part of the Jacobian (only done once, as soon as the kinematics is set)
      virtual void preparePhaseSpace() = 0;
      /// \f$k_{\rm T}\f$-factorised matrix element (event weight)
      /// \return Weight of the point in the phase space to the integral
      virtual double computeKTFactorisedMatrixElement() = 0;
      /// Set the kinematics of the incoming and outgoing protons (or remnants)
      void fillPrimaryParticlesKinematics();
      /// Set the kinematics of the outgoing central system
      virtual void fillCentralParticlesKinematics() = 0;
      /// Set the list of intermediate partons in the process
      void setIntermediatePartons(const std::array<pdgid_t, 2>& part) { intermediate_parts_ = part; }
      /// Set the list of central particles produced
      void setProducedParticles(const std::vector<pdgid_t>& prod) { produced_parts_ = prod; }

      /// Log-virtuality range of the intermediate parton
      Limits log_qt_limits_;
      /// Intermediate azimuthal angle range
      Limits phi_qt_limits_;
      /// Invariant mass range for the scattered excited system
      Limits mx_limits_;

      /// Virtuality of the first intermediate parton (photon, pomeron, ...)
      double qt1_{0.};
      /// Azimuthal rotation of the first intermediate parton's transverse virtuality
      double phi_qt1_{0.};
      /// Virtuality of the second intermediate parton (photon, pomeron, ...)
      double qt2_{0.};
      /// Azimuthal rotation of the second intermediate parton's transverse virtuality
      double phi_qt2_{0.};

      /// First outgoing proton
      Momentum pX_;
      /// Second outgoing proton
      Momentum pY_;

    private:
      /// First and second intermediate parton (photon, pomeron, ...)
      std::array<pdgid_t, 2> intermediate_parts_;
      /// Type of particles produced in the final state
      std::vector<pdgid_t> produced_parts_;
    };
  }  // namespace proc
}  // namespace cepgen

#endif
