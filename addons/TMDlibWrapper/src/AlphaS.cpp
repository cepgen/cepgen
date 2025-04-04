/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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

#include <tmdlib/TMDlib.h>

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"

using namespace cepgen;
using namespace std::string_literals;

class TMDAlphaS final : public Coupling {
public:
  explicit TMDAlphaS(const ParametersList& params) : Coupling(params) {
    tmd_.setVerbosity(steer<int>("verbosity"));
    if (const auto replica = steer<int>("replica"); replica >= 0)
      tmd_.TMDinit(steer<std::string>("name"), replica);
    else
      tmd_.TMDinit(steer<std::string>("name"));
  }

  static ParametersDescription description() {
    auto desc = Coupling::description();
    desc.setDescription("TMDlib alpha(S) evolution algorithm");
    desc.add("name", "PV17_grid_pdf"s).setDescription("dataset name");
    desc.add("verbosity", 99).setDescription("TMDlib evaluator verbosity");
    desc.add("replica", -1).setDescription("dataset replica");
    return desc;
  }

  double operator()(double q) const override { return tmd_.TMDalphas(q); }

private:
  mutable TMDlib::TMD tmd_;
};
REGISTER_ALPHAS_MODULE("tmd", TMDAlphaS);
