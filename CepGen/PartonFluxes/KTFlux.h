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

#include "CepGen/PartonFluxes/PartonFlux.h"

namespace cepgen {
  class KTFlux : public PartonFlux {
  public:
    explicit KTFlux(const ParametersList& params) : PartonFlux(params) {}

    static ParametersDescription description() {
      auto desc = PartonFlux::description();
      desc.setDescription("kT-factorised flux");
      return desc;
    }

    bool ktFactorised() const override final { return true; }

  protected:
    /// Minimal value taken for a \f$\k_{\rm T}\f$-factorised flux
    static constexpr double kMinKTFlux = 1.e-20;
  };
}  // namespace cepgen
