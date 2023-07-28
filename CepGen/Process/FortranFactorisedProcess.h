/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2023  Laurent Forthomme
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

#ifndef CepGen_Process_FortranFactorisedProcess_h
#define CepGen_Process_FortranFactorisedProcess_h

#include <functional>

#include "CepGen/Process/FactorisedProcess.h"

namespace cepgen {
  namespace proc {
    /// Compute the matrix element for a generic collinear, or \f$k_{\rm T}\f$-factorised process defined in a Fortran weighting function
    class FortranFactorisedProcess : public FactorisedProcess {
    public:
      explicit FortranFactorisedProcess(const ParametersList& params, std::function<double(void)> func);
      ProcessPtr clone() const override { return ProcessPtr(new FortranFactorisedProcess(*this)); }

      static ParametersList kProcParameters;

    private:
      void preparePhaseSpace() override;
      double computeFactorisedMatrixElement() override;
      void fillCentralParticlesKinematics() override;

      std::function<double(void)> func_;  ///< Function to be called for weight computation

      // mapped variables
      double m_y1_;           ///< First outgoing particle rapidity
      double m_y2_;           ///< Second outgoing particle rapidity
      double m_pt_diff_;      ///< Transverse momentum balance between outgoing particles
      double m_phi_pt_diff_;  ///< Azimuthal angle difference between outgoing particles
    };
  }  // namespace proc
}  // namespace cepgen

#endif
