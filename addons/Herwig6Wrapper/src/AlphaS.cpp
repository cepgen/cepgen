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

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGenHerwig6/Herwig6Interface.h"

namespace cepgen::herwig6 {
  class AlphaS final : public cepgen::Coupling {
  public:
    explicit AlphaS(const ParametersList& params) : cepgen::Coupling(params), mode_(steer<int>("mode")) {
      if (mode_ < 1 || mode_ > 3)
        throw CG_FATAL("herwig6:AlphaS") << "Invalid mode steered: should be between 1 and 3, got " << mode_ << ".";
      hwpram_.ncolo = steer<int>("ncolo");
      hwpram_.qcdlam = steer<double>("qcdlam");
      hwpram_.qcdl5 = steer<double>("qcdl5");
      hwualf(0, 0.);
    }

    inline static ParametersDescription description() {
      auto desc = cepgen::Coupling::description();
      desc.setDescription("Herwig6 modelling of alpha(S) running");
      initialise();
      desc.add<int>("mode", 1)
          .setDescription("running mode")
          .allow(1, "two-loop flavour thresholds")
          .allow(2, "ratio of mode-1 with 5-flavour beta with Lambda=QCDL3")
          .allow(3, "one-loop with 5-flavour beta and Lambda=QCDL3");
      desc.add<int>("ncolo", hwpram_.ncolo).setDescription("number of colours to consider");
      desc.add<double>("qcdlam", hwpram_.qcdlam).setDescription("5-flavour Lambda_MS-bar at large x/z");
      desc.add<double>("qcdl5", hwpram_.qcdl5).setDescription("5-flavour Lambda_MC");
      return desc;
    }

    inline double operator()(double q) const override { return hwualf(mode_, q * q); }

  private:
    const int mode_;
  };
}  // namespace cepgen::herwig6
using Herwig6AlphaS = cepgen::herwig6::AlphaS;
REGISTER_ALPHAS_MODULE("herwig6", Herwig6AlphaS);
