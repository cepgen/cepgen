/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2025  Laurent Forthomme
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

#include "CepGen/Core/ParametersList.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

EventModifier::EventModifier(const ParametersList& params)
    : EventHandler(params), seed_(steerAs<int, long long>("seed")), max_trials_(steer<int>("maxTrials")) {
  CG_DEBUG("EventModifier:init") << "\"" << name_ << "\"-type event modifier built with:\n\t"
                                 << "* seed = " << seed_ << "\n\t"
                                 << "* maximum trials: " << max_trials_;
}

void EventModifier::readStrings(const std::vector<std::string>& parameters_list) {
  if (parameters_list.empty())
    return;
  for (const auto& p : parameters_list)
    readString(p);
  CG_DEBUG("EventModifier:configure") << "Feeding '" << name_ << "' event modifier algorithm with:\n\t  "
                                      << utils::merge(parameters_list, "\n\t  ");
}

ParametersDescription EventModifier::description() {
  auto desc = EventHandler::description();
  desc.add("seed", -1).setDescription("Random number generator seed");
  desc.add("maxTrials", 1)
      .setDescription("Maximum number of attempts to modify the event before giving up and returning a zero-weight");
  return desc;
}
