/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  const double Kinematics::MX_MIN = 1.07;  // mp+mpi+-

  Kinematics::Kinematics(const ParametersList& params)
      : SteeredObject(params), incoming_beams_(params_), cuts_(params_) {
    CG_DEBUG("Kinematics") << "Building a Kinematics parameters container "
                           << "with the following parameters:\n\t" << params << ".";
    //----- outgoing particles definition
    if (params_.has<std::vector<int> >("minFinalState"))
      for (const auto& pdg : steer<std::vector<int> >("minFinalState"))
        minimum_final_state_.emplace_back((pdgid_t)pdg);

    //--- specify where to look for the grid path for gluon emission
    if (params.has<std::string>("kmrGridPath"))
      kmr::GluonGrid::get(ParametersList(params_).set<std::string>("path", steer<std::string>("kmrGridPath")));
  }

  void Kinematics::setParameters(const ParametersList& params) {
    SteeredObject::setParameters(params);
    incoming_beams_.setParameters(params_);
    cuts_.setParameters(params_);
  }

  ParametersList Kinematics::parameters(bool extended) const {
    ParametersList params;
    params += incoming_beams_.parameters();         // beam particles
    params += cuts_.initial.parameters(extended);   // incoming partons
    params += cuts_.central.parameters(extended);   // central particles
    params += cuts_.remnants.parameters(extended);  // beam remnants
    // minimum final state content
    if (!minimum_final_state_.empty()) {
      std::vector<int> min_pdgs;
      std::transform(
          minimum_final_state_.begin(), minimum_final_state_.end(), std::back_inserter(min_pdgs), [](const auto& pdg) {
            return (int)pdg;
          });
      params.set<std::vector<int> >("minFinalState", min_pdgs);
    }
    // per-PDGid selection
    if (!cuts_.central_particles.empty()) {
      ParametersList per_part;
      for (const auto& cuts_vs_part : cuts_.central_particles)
        per_part.set<ParametersList>(std::to_string(cuts_vs_part.first), cuts_vs_part.second.parameters(extended));
      params.set<ParametersList>("cuts", per_part);
    }
    CG_DEBUG("Kinematics:parameters") << "Kinematics parameters values retrieved:\n"
                                      << ParametersDescription(params) << ".";
    return params;
  }

  ParametersDescription Kinematics::description() {
    auto desc = ParametersDescription();
    desc += IncomingBeams::description();
    desc += CutsList::description();
    return desc;
  }
}  // namespace cepgen
