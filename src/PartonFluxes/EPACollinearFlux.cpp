/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
 *                2009-2012  Nicolas Schul, Jerome de Favereau de Jeneret
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

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/CollinearFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"

using namespace cepgen;

class EPACollinearFlux final : public CollinearFlux {
public:
  explicit EPACollinearFlux(const ParametersList& params)
      : CollinearFlux(params),
        ff_(FormFactorsFactory::get().build(steer<ParametersList>("formFactors"))),
        beam_particle_mass2_(std::pow(steer<ParticleProperties>("pdgId").mass, 2)) {}

  static ParametersDescription description() {
    auto desc = CollinearFlux::description();
    desc.setDescription("EPA FF-dependent flux");
    desc.add("formFactors", FormFactorsFactory::get().describeParameters("StandardDipole"));
    desc.addAs<int, pdgid_t>("pdgId", PDG::proton).setDescription("beam particle PDG id");
    return desc;
  }

  bool fragmenting() const override { return ff_->fragmenting(); }
  pdgid_t partonPdgId() const override { return PDG::photon; }

  double fluxQ2(double x, double q2) const override {
    if (!x_range_.contains(x, true))
      return 0.;
    const auto q2min = utils::kt::q2(x, 0., mass2());
    if (q2min == 0. || q2 < q2min)
      return 0.;
    const auto form_factors = (*ff_)(q2);
    return alpha_over_pi_ * ((1. - x) * (1. - q2min / q2) * form_factors.FE + 0.5 * x * x * form_factors.FM);
  }

private:
  double mass2() const override { return beam_particle_mass2_; }
  const std::unique_ptr<formfac::Parameterisation> ff_;
  const double beam_particle_mass2_;
};
REGISTER_COLLINEAR_FLUX("EPAFlux", EPACollinearFlux);
