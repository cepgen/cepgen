/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2023  Laurent Forthomme
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

#ifndef CepGen_Process_KTPhaseSpaceGenerator_h
#define CepGen_Process_KTPhaseSpaceGenerator_h

#include "CepGen/Process/PhaseSpaceGenerator.h"

namespace cepgen {
  namespace proc {
    /// \f$k_{\rm T}\f$-factorisation phase space generator
    /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
    /// \date Apr 2016
    class KTPhaseSpaceGenerator final : public PhaseSpaceGenerator {
    public:
      explicit KTPhaseSpaceGenerator(Process*);

      bool ktFactorised() const override { return true; }

      void initialise() override;
      bool generatePartonKinematics() override;
      double fluxes() const override;

    protected:
      // mapped variables
      double m_qt1_{0.};      ///< Virtuality of the first intermediate parton (photon, pomeron, ...)
      double m_phi_qt1_{0.};  ///< Azimuthal rotation of the first intermediate parton's transverse virtuality
      double m_qt2_{0.};      ///< Virtuality of the second intermediate parton (photon, pomeron, ...)
      double m_phi_qt2_{0.};  ///< Azimuthal rotation of the second intermediate parton's transverse virtuality
    };
  }  // namespace proc
}  // namespace cepgen

#endif
