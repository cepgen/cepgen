/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2026  Laurent Forthomme
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

#include <array>
#include <memory>

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Process/PartonSampler.h"

namespace cepgen::proc {
  class FactorisedProcess;
}  // namespace cepgen::proc

namespace cepgen {
  class PartonFlux;
  /// Two-parton phase space generator, built from one single-beam parton sampler per incoming beam
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jul 2023
  class PartonsPhaseSpaceGenerator : public SteeredObject<PartonsPhaseSpaceGenerator> {
  public:
    /// Build a two-parton generator from two (positive-z and negative-z) single-beam samplers
    explicit PartonsPhaseSpaceGenerator(const ParametersList&,
                                        std::unique_ptr<PartonSampler> positive_sampler,
                                        std::unique_ptr<PartonSampler> negative_sampler);

    static ParametersDescription description();

    void initialise(proc::FactorisedProcess*);  ///< Build the fluxes and register all integration variables

    bool ktFactorised() const;        ///< Do both incoming partons carry a primordial kT?
    bool generatePartonKinematics();  ///< Generate the 4-momentum of incoming partons
    double fluxes() const;            ///< Retrieve the event weight in the phase space

    /// Retrieve a type-casted positive-z parton flux modelling
    template <typename T = PartonFlux>
    inline const T& positiveFlux() const {
      return sampler(PartonSampler::Side::positive).flux<T>();
    }
    /// Retrieve a type-casted negative-z parton flux modelling
    template <typename T = PartonFlux>
    inline const T& negativeFlux() const {
      return sampler(PartonSampler::Side::negative).flux<T>();
    }
    const PartonSampler& sampler(PartonSampler::Side) const;  ///< Single-beam parton sampler

  private:
    proc::FactorisedProcess& process();              ///< Consumer process object
    const proc::FactorisedProcess& process() const;  ///< Const-qualified consumer process object

    std::array<std::unique_ptr<PartonSampler>, 2> samplers_;  ///< Positive- and negative-z single-beam samplers
    proc::FactorisedProcess* proc_{nullptr};                  //NOT owning
  };
}  // namespace cepgen

#endif
