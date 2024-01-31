/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2023  Laurent Forthomme
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

#ifndef CepGen_Process_Process2to4_h
#define CepGen_Process_Process2to4_h

#include "CepGen/Process/FactorisedProcess.h"

namespace cepgen {
  namespace proc {
    /// A 2-to-4 (or 2-to-2 central) process
    class Process2to4 : public FactorisedProcess {
    public:
      /// Initialise a 2-to-4 process
      /// \param[in] params Collection of user-defined steering parameters
      /// \param[in] cs_id Central particles PDG id
      explicit Process2to4(const ParametersList& params, pdgid_t cs_id);

    protected:
      /// Set all cuts for the single outgoing particle phase space definition
      void setCuts(const cuts::Central& single);

      void prepareFactorisedPhaseSpace() override;
      double computeFactorisedMatrixElement() override;
      void fillCentralParticlesKinematics() override;

      /// Conform all kinematics variables to the user-defined phase space
      virtual void prepareProcessKinematics() = 0;
      /// Computation rule for the central matrix element
      virtual double computeCentralMatrixElement() const = 0;

      //--- Mandelstam variables
      double that() const;  ///< \f$\hat t=\frac{1}{2}\left[(p_1-p_3)^2+(p_2-p_4)^2\right]\f$
      double uhat() const;  ///< \f$\hat u=\frac{1}{2}\left[(p_1-p_4)^2+(p_2-p_3)^2\right]\f$

      ParticleProperties cs_prop_;  ///< PDG properties of the central outgoing particles

      cuts::Central single_limits_;  ///< Limits to be applied on single central system's particles

      // mapped variables
      double m_y_c1_{0.};         ///< Rapidity of the first central particle
      double m_y_c2_{0.};         ///< Rapidity of the second central particle
      double m_pt_diff_{0.};      ///< Transverse momentum difference for the two central particle
      double m_phi_pt_diff_{0.};  ///< Azimuthal angle difference for the two central particles

    private:
      // factor 1/4 from jacobian of transformations
      static constexpr double prefactor_ = 0.25 * 0.0625 * M_1_PI * M_1_PI;
    };
  }  // namespace proc
}  // namespace cepgen

#endif
