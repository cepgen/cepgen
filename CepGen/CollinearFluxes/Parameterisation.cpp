/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2022  Laurent Forthomme
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

#include "CepGen/CollinearFluxes/Parameterisation.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  namespace collflux {
    Parameterisation::Parameterisation(const ParametersList& params)
        : NamedModule(params),
          mp_(PDG::get().mass(PDG::proton)),
          mp2_(mp_ * mp_),
          q2_range_(steer<Limits>("q2range")),
          qscale_(steer<double>("qscale")) {
      q2_range_.min() = std::max(q2_range_.min(), 0.);
    }

    std::string Parameterisation::describe() const { return name_; }

    ParametersDescription Parameterisation::description() {
      auto desc = ParametersDescription();
      desc.setDescription("Unnamed collinear flux");
      desc.add<Limits>("q2range", {0., 1.e4});
      desc.add<double>("qscale", 0.71);
      return desc;
    }

    std::ostream& operator<<(std::ostream& os, const Parameterisation& sf) { return os << sf.describe(); }
  }  // namespace collflux
}  // namespace cepgen
