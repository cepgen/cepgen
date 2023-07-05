/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *                     2003  Maarten Boonekamp, Tibor Kucs
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
        : NamedModule<int>(params),
          mp_(PDG::get().mass(PDG::proton)),
          mp2_(mp_ * mp_),
          t_range_(steer<Limits>("q2range")),
          qscale_(steer<double>("qscale")) {}

    std::string Parameterisation::describe() const {
      std::ostringstream os;
      os << (Type)name_;
      return os.str();
    }

    ParametersDescription Parameterisation::description() {
      auto desc = ParametersDescription();
      desc.setDescription("Unnamed collinear flux");
      desc.add<Limits>("q2range", {0., 1.e4});
      desc.add<double>("qscale", 0.71);
      return desc;
    }

    std::ostream& operator<<(std::ostream& os, const Parameterisation& sf) { return os << sf.describe(); }

    /// Human-readable format of a collinear flux type
    std::ostream& operator<<(std::ostream& os, const collflux::Type& flux) {
      switch (flux) {}
      return os << "<invalid>";
    }
  }  // namespace collflux
}  // namespace cepgen
