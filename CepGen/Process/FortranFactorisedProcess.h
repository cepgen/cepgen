/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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

#include "CepGen/Process/FactorisedProcess.h"

namespace cepgen::proc {
  /// Compute the matrix element for a generic factorised process defined in a Fortran weighting function
  class FortranFactorisedProcess : public FactorisedProcess {
  public:
    /// Construct a Fortran-CepGen interface object using a double precision argument-less F77 function
    /// \param[in] func a double precision argument-less Fortran function returning the event weight
    explicit FortranFactorisedProcess(const ParametersList&, const std::function<double(void)>& func);
    ProcessPtr clone() const override { return ProcessPtr(new FortranFactorisedProcess(*this)); }

    static ParametersList kProcParameters;

  private:
    void prepareFactorisedPhaseSpace() override final;
    double computeFactorisedMatrixElement() override final;

    const std::function<double(void)> func_;  ///< Function to be called for weight computation

    // mapped variables
    double m_y1_{0.};           ///< First outgoing particle rapidity
    double m_y2_{0.};           ///< Second outgoing particle rapidity
    double m_pt_diff_{0.};      ///< Transverse momentum balance between outgoing particles
    double m_phi_pt_diff_{0.};  ///< Azimuthal angle difference between outgoing particles
  };
}  // namespace cepgen::proc

#endif
