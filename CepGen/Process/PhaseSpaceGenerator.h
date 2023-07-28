/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

namespace cepgen {
  class PartonFlux;
  namespace proc {
    class Process;
    /**
     * A generic phase space integration wrapper.
     * \brief Class template to define any phase space helper process
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2023
     */
    class PhaseSpaceGenerator {
    public:
      /// Class constructor
      /// \param[in] params Parameters list
      /// \param[in] output Produced final state particles
      explicit PhaseSpaceGenerator(Process* proc) : proc_(*proc) {}

      const PartonFlux& positiveFlux() const { return *pos_flux_; }
      const PartonFlux& negativeFlux() const { return *neg_flux_; }

      /// Is the phase space generator intended to generate primordial kT for the incoming partons?
      virtual bool ktFactorised() const = 0;

      /// Initialise the process and define the integration phase space
      virtual void init() = 0;
      /// Generate the 4-momentum of incoming partons
      virtual bool generatePartonKinematics() = 0;
      /// Retrieve the event weight in the phase space
      virtual double fluxes() const = 0;

      /// Consumer process object
      Process& process() { return proc_; }
      /// Const-qualified consumer process object
      const Process& process() const { return const_cast<const Process&>(proc_); }

    protected:
      std::unique_ptr<PartonFlux> pos_flux_{nullptr}, neg_flux_{nullptr};

    private:
      Process& proc_;  //NOT owning
    };
  }  // namespace proc
}  // namespace cepgen

#endif
