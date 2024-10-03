/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Message.h"
#include "CepGenHerwig6/Herwig6Interface.h"

namespace cepgen::herwig6 {
  class StructureFunctions final : public strfun::Parameterisation {
  public:
    explicit StructureFunctions(const ParametersList& params)
        : strfun::Parameterisation(params), idhad_(steer<int>("idhad")), nset_(steer<int>("nset")) {}

    inline static ParametersDescription description() {
      initialise();
      auto desc = strfun::Parameterisation::description();
      desc.setDescription("Herwig 6 structure functions evaluator");
      desc.add<int>("idhad", 73)
          .setDescription("type of hadron")
          .allow(30, "pi-")
          .allow(38, "pi+")
          .allow(59, "photon")
          .allow(73, "proton")
          .allow(75, "neutron")
          .allow(91, "antiproton")
          .allow(93, "antineutron");
      desc.add<int>("nset", 8)
          .setDescription("structure functions set")
          .allow(1, "Duke & Owens set 1 (for soft/hard glue)")
          .allow(2, "Duke & Owens set 2 (for soft/hard glue)")
          .allow(3, "Eichten & al. set 1 (nucleons only)")
          .allow(4, "Eichten & al. set 2 (nucleons only)")
          .allow(5, "Owens set 1.1")
          .allow(6, "MRST98LO (central alpha(S)/gluon)")
          .allow(7, "MRST98LO (higher gluon)")
          .allow(8, "MRST98LO (average of central and higher gluon)");
      return desc;
    }

  private:
    inline void eval() override { setF2(hwsfun(args_.xbj, args_.q2, idhad_, nset_, 2)); }
    const int idhad_, nset_;
  };
}  // namespace cepgen::herwig6
using Herwig6StructureFunctions = cepgen::herwig6::StructureFunctions;
REGISTER_STRFUN("herwig6", 403, Herwig6StructureFunctions);
