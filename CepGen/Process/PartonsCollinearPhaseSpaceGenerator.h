/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

#ifndef CepGen_Process_PartonsCollinearPhaseSpaceGenerator_h
#define CepGen_Process_PartonsCollinearPhaseSpaceGenerator_h

#include "CepGen/Process/PartonsPhaseSpaceGenerator.h"

namespace cepgen {
  namespace proc {
    /// Collinear factorisation phase space generator
    /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
    /// \date Jul 2023
    class PartonsCollinearPhaseSpaceGenerator final : public PartonsPhaseSpaceGenerator {
    public:
      explicit PartonsCollinearPhaseSpaceGenerator(FactorisedProcess*);

      bool ktFactorised() const override { return false; }

      void initialise() override;
      bool generatePartonKinematics() override;
      double fluxes() const override;

    protected:
      // mapped variables
      double m_t1_{0.}, m_t2_{0.};
    };
  }  // namespace proc
}  // namespace cepgen

#endif
