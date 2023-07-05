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

    double operator()(double, double = 0.) const override { return 0.; }
    bool ktFactorised() const override final { return true; }

  protected:
    /// Minimal value taken for a \f$\k_{\rm T}\f$-factorised flux
    static constexpr double kMinKTFlux = 1.e-20;
    virtual double mass2() const { return 0.; }

    struct Q2Values {
      double min{0.}, q2{0.};
    };
    /// Compute the minimum and kT-dependent Q^2
    Q2Values computeQ2(double x, double kt2, double mx2 = 0.) const {
      Q2Values out;
      const auto mi2 = mass2();
      const auto dm2 = (mx2 == 0.) ? 0. : mx2 - mi2;
      out.min = ((x * dm2) + x * x * mi2) / (1. - x);
      out.q2 = out.min + kt2 / (1. - x);
      return out;
    }
  };
}  // namespace cepgen

