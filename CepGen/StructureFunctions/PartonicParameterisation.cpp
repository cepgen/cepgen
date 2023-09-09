/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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
#include "CepGen/StructureFunctions/PartonicParameterisation.h"

namespace cepgen {
  namespace strfun {
    constexpr std::array<short, 6> PartonicParameterisation::Q_TIMES_3, PartonicParameterisation::QUARK_PDGS;

    PartonicParameterisation::PartonicParameterisation(const ParametersList& params)
        : Parameterisation(params), num_flavours_(steer<int>("numFlavours")), mode_(steerAs<int, Mode>("mode")) {
      if (num_flavours_ == 0 || num_flavours_ > QUARK_PDGS.size())
        throw CG_FATAL("Partonic") << "Invalid number of flavours (" << num_flavours_ << " selected.";
    }

    ParametersDescription PartonicParameterisation::description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Partonic structure functions parameterisation");
      desc.add<int>("numFlavours", 4).setDescription("Number of parton flavours to consider in summation");
      desc.add<int>("mode", (int)Mode::full);
      return desc;
    }

    void PartonicParameterisation::eval() {
      double f2 = 0.;
      for (int i = 0; i < num_flavours_; ++i) {
        const double prefactor = 1. / 9. * Q_TIMES_3.at(i) * Q_TIMES_3.at(i);
        const double xq = evalxQ2(QUARK_PDGS.at(i), args_.xbj, args_.q2),
                     xqbar = evalxQ2(-QUARK_PDGS.at(i), args_.xbj, args_.q2);
        switch (mode_) {
          case Mode::full:
            f2 += prefactor * (xq + xqbar);
            break;
          case Mode::valence:
            f2 += prefactor * (xq - xqbar);
            break;
          case Mode::sea:
            f2 += prefactor * (2. * xqbar);
            break;
        }
      }
      setF2(f2);
    }

    std::ostream& operator<<(std::ostream& os, const strfun::PartonicParameterisation::Mode& mode) {
      switch (mode) {
        case strfun::PartonicParameterisation::Mode::full:
          return os << "all quarks";
        case strfun::PartonicParameterisation::Mode::valence:
          return os << "valence quarks";
        case strfun::PartonicParameterisation::Mode::sea:
          return os << "sea quarks";
      }
      return os;
    }
  }  // namespace strfun
}  // namespace cepgen
