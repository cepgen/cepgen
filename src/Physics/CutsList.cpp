/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2025  Laurent Forthomme
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

#include "CepGen/Physics/CutsList.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

CutsList::CutsList(const ParametersList& params)
    : SteeredObject(params),
      initial(steer<ParametersList>("initial")),
      central(steer<ParametersList>("central")),
      remnants(steer<ParametersList>("remnants")) {}

void CutsList::setParameters(const ParametersList& params) {
  if (params.empty())
    return;
  initial.setDescribedParameters(params);
  central.setDescribedParameters(params);
  remnants.setDescribedParameters(params);
  if (params.has<ParametersList>("cuts")) {  // per-particle cuts
    const auto per_part_cuts = params.get<ParametersList>("cuts");
    for (const auto& part : per_part_cuts.keys())
      central_particles[static_cast<pdgid_t>(std::stoi(part))].setDescribedParameters(
          per_part_cuts.get<ParametersList>(part));
  }

  // override the parameters from sub-parameters content
  params_.set("initial", initial.parameters())
      .set("central", central.parameters())
      .set("remnants", remnants.parameters());
  for (const auto& cuts_vs_part : central_particles)  // per-PDGid selection
    params_.operator[]<ParametersList>("cuts").set(std::to_string(cuts_vs_part.first),
                                                   cuts_vs_part.second.parameters());
  CG_DEBUG("CutsList:setParameters") << "User specified the following cuts list:\n" << *this << ".";
  for (const auto& key : params_.keys()) {
    if (key == "initial" || key == "central" || key == "remnants" || key == "cuts") {
      auto& cuts = params_.operator[]<ParametersList>(key);
      for (const auto& lim_key : cuts.keysOf<Limits>())
        if (auto& lim = cuts.operator[]<Limits>(lim_key); lim.min() == 0. && lim.max() == 0.) {
          CG_WARNING("CutsList:setParameters")
              << "Unset the range for '" << key << "/" << lim_key << "' from " << lim << ".";
          lim = Limits{};
        }
    } else
      params_.erase(key);
  }
  initial.setParameters(steer<ParametersList>("initial"));
  central.setParameters(steer<ParametersList>("central"));
  remnants.setParameters(steer<ParametersList>("remnants"));
}

ParametersDescription CutsList::description() {
  auto desc = ParametersDescription();
  desc.add("initial", cuts::Initial::description());
  desc.add("central", cuts::Central::description());
  desc.add("remnants", cuts::Remnants::description());
  return desc;
}

namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const CutsList& cl) {
    auto dump_cuts = [&os](const auto& obj) {
      std::string sep;
      for (const auto& lim : obj.parameters().template keysOf<Limits>()) {
        const auto& limit = obj.parameters().template get<Limits>(lim);
        if (limit.valid() && obj.description().has(lim))
          os << sep << obj.description().get(lim).description() << ": " << limit, sep = ";";
      }
    };
    os << "init.system{";
    dump_cuts(cl.initial);
    os << "}, cent.system{";
    dump_cuts(cl.central);
    os << "}, remnants{";
    dump_cuts(cl.remnants);
    return os << "}";
  }
}  // namespace cepgen
