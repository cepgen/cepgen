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

#include "CepGen/Modules/NamedModule.h"

namespace cepgen::proc {
  class FactorisedProcess;
}
namespace cepgen::cuts {
  struct Central;
}

namespace cepgen {
  /// Class template to define any phase space helper process
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Feb 2024
  class PhaseSpaceGenerator : public NamedModule<PhaseSpaceGenerator> {
  public:
    explicit PhaseSpaceGenerator(const ParametersList& params) : NamedModule(params) {}

    inline virtual bool ktFactorised() const { return false; }

    inline virtual void setCentralCuts(const cuts::Central&) {}  ///< Set cuts on central particles
    virtual void initialise(proc::FactorisedProcess*) = 0;       ///< Set all process parameters

    virtual bool generate() = 0;        ///< Generate a kinematics combination, and return a success flag
    virtual double weight() const = 0;  ///< Return the event weight for a kinematics combination

    virtual pdgids_t partons() const = 0;                  ///< List of incoming partons in kinematics
    virtual void setCentral(const std::vector<int>&) = 0;  ///< Override the central particles list
    virtual std::vector<int> central() const = 0;          ///< List of outgoing central particles in kinematics

    // Mandelstam variables
    virtual double that() const = 0;
    virtual double uhat() const = 0;
  };
}  // namespace cepgen

#endif
