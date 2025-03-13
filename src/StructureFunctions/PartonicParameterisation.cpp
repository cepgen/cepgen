/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

using namespace cepgen;
using namespace cepgen::strfun;

constexpr std::array<short, 6> PartonicParameterisation::Q_TIMES_3, PartonicParameterisation::QUARK_PDG_IDS;

PartonicParameterisation::PartonicParameterisation(const ParametersList& params)
    : Parameterisation(params), num_flavours_(steer<int>("numFlavours")), mode_(steerAs<int, Mode>("mode")) {
  if (num_flavours_ == 0 || num_flavours_ > QUARK_PDG_IDS.size())
    throw CG_FATAL("Partonic") << "Invalid number of flavours (" << num_flavours_ << " selected.";
}

ParametersDescription PartonicParameterisation::description() {
  auto desc = Parameterisation::description();
  desc.setDescription("Partonic structure functions parameterisation");
  desc.add<int>("numFlavours", 4)
      .setDescription("Number of parton flavours to consider in summation")
      .allow(1, "down quark only")
      .allow(2, "down+up quarks")
      .allow(3, "down+up+strange quarks")
      .allow(4, "down+up+strange+charm quarks")
      .allow(5, "down+up+strange+charm+bottom quarks")
      .allow(6, "all quark flavours");
  desc.addAs<int, Mode>("mode", Mode::full);
  return desc;
}

void PartonicParameterisation::eval() {
  double f2 = 0.;
  for (int i = 0; i < num_flavours_; ++i) {
    const double prefactor = 1. / 9. * Q_TIMES_3.at(i) * Q_TIMES_3.at(i);
    const double xq = evalxQ2(QUARK_PDG_IDS.at(i), args_.xbj, args_.q2),
                 xqbar = evalxQ2(-QUARK_PDG_IDS.at(i), args_.xbj, args_.q2);
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

namespace cepgen::strfun {
  std::ostream& operator<<(std::ostream& os, const PartonicParameterisation::Mode& mode) {
    switch (mode) {
      case PartonicParameterisation::Mode::full:
        return os << "all quarks";
      case PartonicParameterisation::Mode::valence:
        return os << "valence quarks";
      case PartonicParameterisation::Mode::sea:
        return os << "sea quarks";
    }
    return os;
  }
}  // namespace cepgen::strfun
