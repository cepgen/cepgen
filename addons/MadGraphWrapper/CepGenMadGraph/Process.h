/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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

#ifndef CepGenMadGraph_Process_h
#define CepGenMadGraph_Process_h

#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Physics/ParticleProperties.h"

namespace cepgen {
  class Momentum;
}  // namespace cepgen

namespace cepgen::mg5amc {
  /// Wrapper around a generic MadGraph process definition
  class Process : public NamedModule<Process> {
  public:
    explicit Process(const ParametersList&);
    virtual ~Process() = default;

    static ParametersDescription description();

    inline const spdgids_t& intermediatePartons() const { return incoming_pdgids_; }
    inline const spdgids_t& centralSystem() const { return central_pdgids_; }

    virtual void initialise(const std::string&) = 0;
    virtual double eval() = 0;
    virtual const std::vector<Momentum>& momenta() = 0;

    Process& setMomentum(size_t i, const Momentum& mom);

  protected:
    const spdgids_t incoming_pdgids_;  ///< incoming partons content
    const spdgids_t central_pdgids_;   ///< central system particles content
    std::vector<double*> mom_;
  };
}  // namespace cepgen::mg5amc

#endif
