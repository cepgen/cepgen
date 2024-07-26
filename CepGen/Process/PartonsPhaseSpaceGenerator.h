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

#ifndef CepGen_Process_PartonsPhaseSpaceGenerator_h
#define CepGen_Process_PartonsPhaseSpaceGenerator_h

#include <memory>

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  class PartonFlux;
  namespace proc {
    class FactorisedProcess;
  }
  /**
     * A generic phase space integration wrapper.
     * \brief Class template to define any phase space helper process
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2023
     */
  class PartonsPhaseSpaceGenerator : public SteeredObject<PartonsPhaseSpaceGenerator> {
  public:
    explicit PartonsPhaseSpaceGenerator(const ParametersList& params) : SteeredObject(params) {}

    void initialise(proc::FactorisedProcess* process) {
      proc_ = process;
      initialise();
    }

    virtual bool ktFactorised() const = 0;        ///< Do incoming partons carry a primordial kT?
    virtual bool generatePartonKinematics() = 0;  ///< Generate the 4-momentum of incoming partons
    virtual double fluxes() const = 0;            ///< Retrieve the event weight in the phase space

    /// Retrieve a type-casted positive-z parton flux modelling
    template <typename T = PartonFlux>
    inline const T& positiveFlux() const {
      return dynamic_cast<const T&>(*pos_flux_);
    }
    /// Retrieve a type-casted negative-z parton flux modelling
    template <typename T = PartonFlux>
    inline const T& negativeFlux() const {
      return dynamic_cast<const T&>(*neg_flux_);
    }

  protected:
    virtual void initialise() = 0;  ///< Initialise the process and define the integration phase space
    inline proc::FactorisedProcess& process() { return *proc_; }  ///< Consumer process object
    /// Const-qualified consumer process object
    inline const proc::FactorisedProcess& process() const { return const_cast<const proc::FactorisedProcess&>(*proc_); }
    std::unique_ptr<PartonFlux> pos_flux_{nullptr}, neg_flux_{nullptr};

  private:
    proc::FactorisedProcess* proc_{nullptr};  //NOT owning
  };
}  // namespace cepgen

#endif
