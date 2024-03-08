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

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/ModuleFactory.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"

namespace cepgen {
  template class ModuleFactory<MadGraphProcess, std::string>;

  MadGraphProcess::MadGraphProcess(const ParametersList& params)
      : SteeredObject(params),
        incoming_pdgids_(steer<std::vector<int> >("incomingSystem").begin(),
                         steer<std::vector<int> >("incomingSystem").end()),
        central_pdgids_(steer<std::vector<int> >("outgoingSystem").begin(),
                        steer<std::vector<int> >("outgoingSystem").end()) {}

  MadGraphProcess& MadGraphProcess::setMomentum(size_t i, const Momentum& mom) {
    if (i > mom_.size())
      throw CG_FATAL("MadGraphProcess") << "Invalid index for momentum: " << i << "!";
    mom_[i][0] = mom.energy();
    mom_[i][1] = mom.px();
    mom_[i][2] = mom.py();
    mom_[i][3] = mom.pz();
    return *this;
  }

  ParametersDescription MadGraphProcess::description() {
    auto desc = ParametersDescription();
    desc.setDescription("generic mg5_aMC@NLO process");
    desc.add<std::vector<int> >("incomingSystem", {}).setDescription("list of incoming partons for the process");
    desc.add<std::vector<int> >("outgoingSystem", {}).setDescription("list of central particles generated");
    return desc;
  }
}  // namespace cepgen
