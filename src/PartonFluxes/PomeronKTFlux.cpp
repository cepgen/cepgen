/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"

using namespace cepgen;

/// Elastic, kT-dependent pomeron emission from a proton (Donnachie-Landshoff modelling)
class PomeronKTFlux final : public KTFlux {
public:
  explicit PomeronKTFlux(const ParametersList& params)
      : KTFlux(params),
        beta2_(steer<double>("beta2")),
        epsilon_(steer<double>("epsilon")),
        alpha_prime_(steer<double>("alphaPrime")),
        f1_scale_(steer<double>("formFactorScale")) {}

  static ParametersDescription description() {
    auto desc = KTFlux::description();
    desc.setDescription("Elastic pomeron emission from a proton (Donnachie-Landshoff)");
    desc.add("beta2", 3.5).setDescription("squared pomeron-quark coupling (in GeV^-2)");
    desc.add("epsilon", 0.0808).setDescription("pomeron trajectory intercept, alpha(0)-1");
    desc.add("alphaPrime", 0.25).setDescription("pomeron trajectory slope (in GeV^-2)");
    desc.add("formFactorScale", 0.71).setDescription("proton isoscalar form factor dipole scale (in GeV^2)");
    return desc;
  }

  bool fragmenting() const override { return false; }
  double mass2() const override { return mp2_; }
  spdgid_t partonPdgId() const override { return PDG::pomeron; }

  double fluxMX2(double x, double kt2, double) const override {
    if (!x_range_.contains(x, true))
      return 0.;
    const auto t = -utils::kt::q2(x, kt2, mass2());  // pomeron four-momentum transfer, t < 0
    const auto alpha = 1. + epsilon_ + alpha_prime_ * t;
    const auto f1 = diracFormFactor(t);
    return 9. * beta2_ * M_1_PI * M_1_PI * 0.25 * f1 * f1 * std::pow(x, 2. * (1. - alpha)) * kt2;
  }

private:
  /// Donnachie-Landshoff isoscalar proton Dirac form factor
  double diracFormFactor(double t) const {
    const auto four_m2 = 4. * mp2_;
    const auto dipole = 1. / (1. - t / f1_scale_);
    return (four_m2 - 2.8 * t) / (four_m2 - t) * dipole * dipole;
  }

  const double beta2_, epsilon_, alpha_prime_, f1_scale_;
};
REGISTER_KT_FLUX("PomeronElastic", 30, PomeronKTFlux);
