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

#include <cmath>

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen::strfun {
  /// F2 parameterisation for Q^2 < 10 GeV^2
  /// \cite Donnachie:1993it
  class DonnachieLandshoff final : public Parameterisation {
  public:
    explicit DonnachieLandshoff(const ParametersList& params)
        : Parameterisation(params),
          pA_(steer<double>("A")),
          pB_(steer<double>("B")),
          pa_(steer<double>("a")),
          pb_(steer<double>("b")),
          epsilon_(steer<double>("epsilon")),
          delta_r_(steer<double>("deltaR")) {}

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Donnachie-Landshoff");
      desc.add<double>("A", 0.324);
      desc.add<double>("B", 0.098);
      desc.add<double>("a", 0.561991692786383);
      desc.add<double>("b", 0.011133);
      desc.add<double>("epsilon", 0.0808);
      desc.add<double>("deltaR", 0.5475);
      return desc;
    }

    inline void eval() override {
      setF2(pA_ * std::pow(args_.xbj, -epsilon_) * std::pow(args_.q2 / (args_.q2 + pa_), 1 + epsilon_) +
            pB_ * std::pow(args_.xbj, 1. - delta_r_) * std::pow(args_.q2 / (args_.q2 + pb_), delta_r_));
    }

  private:
    const double pA_, pB_, pa_, pb_, epsilon_, delta_r_;
  };
}  // namespace cepgen::strfun
using cepgen::strfun::DonnachieLandshoff;
REGISTER_STRFUN("DonnachieLandshoff", 105, DonnachieLandshoff);
