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

#ifndef CepGen_Process_CentralPhaseSpaceGenerator_h
#define CepGen_Process_CentralPhaseSpaceGenerator_h

#include <memory>

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Physics/Cuts.h"

namespace cepgen {
  namespace proc {
    class FactorisedProcess;
  }
  /// Generic central kinematics generator
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Feb 2024
  class CentralPhaseSpaceGenerator : public SteeredObject<CentralPhaseSpaceGenerator> {
  public:
    explicit CentralPhaseSpaceGenerator(const ParametersList& params) : SteeredObject(params) {}

    void initialise(proc::FactorisedProcess* process) {
      proc_ = process;
      initialise();
    }
    /// Number of dimensions required to generate the kinematics
    virtual size_t ndim() const = 0;
    /// List of produced particles PDG id
    virtual const pdgids_t& particles() const = 0;

    /// Generate the 4-momentum of incoming partons
    virtual double generateKinematics() = 0;

    /// Set all cuts for the single outgoing particle phase space definition
    void setCuts(const cuts::Central& single) { single_limits_ = single; }

  protected:
    static constexpr double NUM_LIMITS = 1.e-3;  ///< Numerical limits for sanity comparisons (MeV/mm-level)
    /// Initialise the process and define the integration phase space
    virtual void initialise() = 0;

    inline proc::FactorisedProcess& process() { return *proc_; }  ///< Consumer process object
    /// Const-qualified consumer process object
    inline const proc::FactorisedProcess& process() const { return const_cast<const proc::FactorisedProcess&>(*proc_); }

    cuts::Central single_limits_;  ///< Limits to be applied on single central system's particles

  private:
    proc::FactorisedProcess* proc_{nullptr};  //NOT owning
  };
}  // namespace cepgen

#endif
