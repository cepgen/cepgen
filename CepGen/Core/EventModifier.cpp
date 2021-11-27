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

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

namespace cepgen {
  EventModifier::EventModifier(const ParametersList& plist)
      : NamedModule(plist),
        seed_(plist.getAs<int, long long>("seed", -1)),
        max_trials_(plist.get<int>("maxTrials", 1)) {
    CG_DEBUG("EventModifier:init") << "\"" << name_ << "\"-type event modifier built with:\n\t"
                                   << "* seed = " << seed_ << "\n\t"
                                   << "* maximum trials: " << max_trials_;
  }

  void EventModifier::readStrings(const std::vector<std::string>& params) {
    if (params.empty())
      return;
    std::ostringstream os;
    for (const auto& p : params) {
      readString(p);
      os << "\n\t  '" << p << "'";
    }
    CG_DEBUG("EventModifier:configure") << "Feeding \"" << name_ << "\" event modifier algorithm with:" << os.str();
  }
}  // namespace cepgen
