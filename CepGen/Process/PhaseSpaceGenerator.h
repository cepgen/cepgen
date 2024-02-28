/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#ifndef CepGen_Process_PhaseSpaceGenerator_h
#define CepGen_Process_PhaseSpaceGenerator_h

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  namespace proc {
    class FactorisedProcess;
  }
  namespace cuts {
    class Central;
  }
  /// Class template to define any phase space helper process
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Feb 2024
  class PhaseSpaceGenerator : public SteeredObject<PhaseSpaceGenerator> {
  public:
    explicit PhaseSpaceGenerator(const ParametersList& params) : SteeredObject(params) {}

    virtual bool ktFactorised() const { return false; }

    virtual void setCentralCuts(const cuts::Central&) const {}  ///< Set cuts on central particles
    virtual void initialise(proc::FactorisedProcess*) = 0;      ///< Set all process parameters
    virtual double generate() = 0;         ///< Generate a kinematics combination and return a weight
    virtual pdgids_t partons() const = 0;  ///< List of incoming partons in kinematics
    virtual pdgids_t central() const = 0;  ///< List of outgoing central particles in kinematics
  };
}  // namespace cepgen

#endif
