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

#include <cmath>

#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/CollinearFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"

using namespace cepgen;

/// Virtuality-dependent Drees-Zeppenfeld photon flux
/// \note Corresponds to PDF:Proton2gammaSet=2 in Pythia 8
/// \cite Drees:1988pp
class DreesZeppenfeldCollinearFlux final : public CollinearFlux {
public:
  explicit DreesZeppenfeldCollinearFlux(const ParametersList& params)
      : CollinearFlux(params),
        scale_(steer<double>("scale")),
        coefficients_a_(steer<std::vector<double> >("coeffsA")) {}

  static ParametersDescription description() {
    auto desc = CollinearFlux::description();
    desc.setDescription("Drees-Zeppenfeld Q^{2}-dependent flux");
    desc.add<double>("scale", 0.71).setDescription("factorisation scale (in GeV^2)");
    desc.add<std::vector<double> >("coeffsA", {-11. / 6, 3., -1.5, 1. / 3});
    return desc;
  }

  bool fragmenting() const override { return true; }
  pdgid_t partonPdgId() const override { return PDG::photon; }

  double fluxQ2(double x, double q2) const override {
    if (!x_range_.contains(x, true))
      return 0.;
    const auto q2min = utils::kt::q2(x, 0. /* kt2 */, mp2_);
    const auto fq4 = std::pow(1 + q2 / scale_, -4);  // Q^2-dependent form factor
    return alpha_over_pi_ * 0.5 * (1. + std::pow(1. - x, 2)) * factorA(1. + scale_ / q2min) * fq4;
  }

private:
  double mass2() const override { return mp2_; }
  double factorA(double a) const {
    auto ret = std::log(a);
    for (int i = 0; i < static_cast<int>(coefficients_a_.size()); ++i)
      ret += coefficients_a_.at(i) * std::pow(a, -i);
    return ret;
  }

  const double scale_;
  const std::vector<double> coefficients_a_;
};
REGISTER_COLLINEAR_FLUX("DreesZeppenfeld", DreesZeppenfeldCollinearFlux);
