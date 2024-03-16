/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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
#include "CepGen/Utils/RandomGenerator.h"

namespace cepgen {
  namespace utils {
    RandomGenerator::RandomGenerator(const ParametersList& params)
        : NamedModule(params), seed_(steer<unsigned long long>("seed")) {}

    double RandomGenerator::exponential(double /*exponent*/) {
      CG_WARNING("RandomGenerator:exponential")
          << "Exponential distribution not implemented for this random number generator.";
      return 0.;
    }

    double RandomGenerator::breitWigner(double /*mean*/, double /*scale*/) {
      CG_WARNING("RandomGenerator:breitWigner")
          << "Breit-Wigner/Cauchy distribution not implemented for this random number generator.";
      return 0.;
    }

    double RandomGenerator::landau(double /*location*/, double /*width*/) {
      CG_WARNING("RandomGenerator:landau") << "Landau distribution not implemented for this random number generator.";
      return 0.;
    }

    int RandomGenerator::poisson(double /*mean*/) {
      CG_WARNING("RandomGenerator:poisson") << "Poisson distribution not implemented for this random number generator.";
      return 0;
    }

    void* RandomGenerator::enginePtr() {
      throw CG_FATAL("RandomGenerator:enginePtr") << "No engine object declared for this random generator.";
    }

    ParametersDescription RandomGenerator::description() {
      auto desc = ParametersDescription();
      desc.setDescription("unnamed random generator");
      desc.add<unsigned long long>("seed", time(nullptr)).setDescription("Random number generator seed");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen
