/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2024  Laurent Forthomme
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

#ifndef CepGenAddOns_MadGraphWrapper_MadGraphProcess_h
#define CepGenAddOns_MadGraphWrapper_MadGraphProcess_h

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/ParticleProperties.h"

namespace cepgen {
  /// Wrapper around a generic MadGraph process definition
  class MadGraphProcess : public SteeredObject<MadGraphProcess> {
  public:
    explicit MadGraphProcess(const ParametersList&);
    virtual ~MadGraphProcess() = default;

    static ParametersDescription description();

    inline const std::vector<int>& intermediatePartons() const { return incoming_pdgids_; }
    inline const std::vector<int>& centralSystem() const { return central_pdgids_; }

    virtual void initialise(const std::string&) = 0;
    virtual double eval() = 0;
    virtual const std::vector<Momentum>& momenta() = 0;

    MadGraphProcess& setMomentum(size_t i, const Momentum& mom);

  protected:
    const std::vector<int> incoming_pdgids_, central_pdgids_;
    std::vector<double*> mom_;
  };
}  // namespace cepgen

#endif
