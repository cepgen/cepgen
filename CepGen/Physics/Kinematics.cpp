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

  Kinematics::Kinematics(const ParametersList& params) : incoming_beams_(params) {
    CG_DEBUG("Kinematics") << "Building a Kinematics parameters container "
                           << "with the following parameters:\n\t" << params << ".";
    //----- phase space definition
    setParameters(params);
  }

  void Kinematics::setParameters(const ParametersList& params) {
    //--- initial partons
    cuts_.initial.setParameters(params);

    //--- central system
    cuts_.central.setParameters(params);
    if (params.has<Limits>("phiptdiff")) {
      CG_WARNING("Kinematics") << "\"phiptdiff\" parameter is deprecated! "
                               << "Please use \"phidiff\" instead.";
      params.fill<Limits>("phiptdiff", cuts_.central.phi_diff());  //legacy
    }
    if (params.has<std::vector<int> >("minFinalState"))
      for (const auto& pdg : params.get<std::vector<int> >("minFinalState"))
        minimum_final_state_.emplace_back((pdgid_t)pdg);
    if (params.has<ParametersList>("cuts")) {  // per-particle cuts
      const auto& per_parts = params.get<ParametersList>("cuts");
      for (const auto& part : per_parts.keys())
        cuts_.central_particles[(pdgid_t)stoi(part)].setParameters(per_parts.get<ParametersList>(part));
    }

    //--- outgoing remnants
    cuts_.remnants.setParameters(params);
    // sanity check
    if (cuts_.remnants.mx().min() < MX_MIN) {
      CG_WARNING("Kinematics:setParameters") << "Minimum diffractive mass set to " << MX_MIN << " GeV.";
      cuts_.remnants.mx().min() = MX_MIN;
    }

    //--- specify where to look for the grid path for gluon emission
    if (params.has<std::string>("kmrGridPath"))
      kmr::GluonGrid::get(params.get<std::string>("kmrGridPath"));
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
    return params;
  }

  std::ostream& operator<<(std::ostream& os, const Kinematics::CutsList& kin) {
    std::string sep;
    os << "initial: {";
    for (const auto& cut : kin.initial.list())
      os << sep << cut, sep = ", ";
    os << "}, central: {";
    sep.clear();
    for (const auto& cut : kin.central.list())
      os << sep << cut, sep = ", ";
    os << "}, remnants: {";
    sep.clear();
    for (const auto& cut : kin.remnants.list())
      os << sep << cut, sep = ", ";
    return os << "}";
  }

  //--------------------------------------------------------------------
  // List of kinematics limits
  //--------------------------------------------------------------------

  Kinematics::CutsList::CutsList()
      : initial(ParametersList().set<Limits>("q2", {0., 1.e5})),
        central(ParametersList().set<double>("ptmin", 0.)),
        remnants(ParametersList().set<Limits>("mx", {MX_MIN, 1000.})) {}
}  // namespace cepgen
