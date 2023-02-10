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

#include <REvolver.h>

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"

namespace cepgen {
  class AlphaSREvolver final : public Coupling {
  public:
    explicit AlphaSREvolver(const ParametersList& params)
        : Coupling(params),
          qc_(steer<double>("qCentral")),
          qevol_(steer<double>("qEvol")),
          order_(steer<int>("order")),
          central2_(revo::RunPar(order_, qc_, qevol_), order_) {}

    static ParametersDescription description() {
      auto desc = Coupling::description();
      desc.setDescription("REvolver alpha(S) evolution algorithm");
      desc.add<int>("order", 5);
      desc.add<double>("qCentral", 0.0822);
      desc.add<double>("qEvol", 1508.04);
      return desc;
    }

    double operator()(double q) const override { return central2_.alpha(q); }

  private:
    const double qc_, qevol_;
    const int order_;
    revo::Core central2_;
  };
}  // namespace cepgen

REGISTER_ALPHAS_MODULE("revolver", AlphaSREvolver);
