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

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/Hadroniser.h"

namespace cepgen {
  namespace hadr {
    Hadroniser::Hadroniser(const ParametersList& plist)
        : EventModifier(plist), remn_fragm_(steer<bool>("remnantsFragmentation")) {}

    ParametersDescription Hadroniser::description() {
      auto desc = EventModifier::description();
      desc.add<bool>("remnantsFragmentation", true)
          .setDescription("Apply the fragmentation algorithm to proton remnants");
      return desc;
    }
  }  // namespace hadr
}  // namespace cepgen
