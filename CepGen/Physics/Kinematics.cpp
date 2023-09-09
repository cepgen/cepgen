/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  Kinematics::Kinematics(const ParametersList& params) : SteeredObject(params) {
    if (!params.empty())
      setParameters(params);
  }

  void Kinematics::setParameters(const ParametersList& params) {
    SteeredObject::setParameters(params);
    CG_DEBUG("Kinematics") << "Building a Kinematics parameters container "
                           << "with the following parameters:\n\t" << params_ << ".";

    incoming_beams_.setParameters(params_);
    cuts_.setParameters(params_);
    //----- outgoing particles definition
    if (params_.has<std::vector<int> >("minFinalState"))
      for (const auto& pdg : steer<std::vector<int> >("minFinalState"))
        minimum_final_state_.emplace_back((pdgid_t)pdg);

    //--- specify where to look for the grid path for gluon emission
    const auto& kmr_grid_path = steerPath("kmrGridPath");
    if (!kmr_grid_path.empty())
      kmr::GluonGrid::get(ParametersList(params_).set<std::string>("path", kmr_grid_path));
  }

  ParametersList Kinematics::parameters(bool extended) const {
    ParametersList params;
    params += incoming_beams_.parameters();  // beam particles
    params += cuts_.parameters(extended);
    // minimum final state content
    if (!minimum_final_state_.empty()) {
      std::vector<int> min_pdgs;
      std::transform(
          minimum_final_state_.begin(), minimum_final_state_.end(), std::back_inserter(min_pdgs), [](const auto& pdg) {
            return (int)pdg;
          });
      params.set<std::vector<int> >("minFinalState", min_pdgs);
    }
    return params;
  }

  ParametersDescription Kinematics::description() {
    auto desc = ParametersDescription();
    desc += IncomingBeams::description();
    desc += CutsList::description();
    desc.add<std::string>("kmrGridPath", "").setDescription("path to the KMR interpolation grid");
    return desc;
  }
}  // namespace cepgen
