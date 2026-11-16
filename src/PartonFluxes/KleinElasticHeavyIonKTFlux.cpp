/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2026  Laurent Forthomme
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
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"

using namespace cepgen;

/// Realistic nuclear form-factor as used in STARLIGHT
/// See \cite Klein:2016yzr
class KleinElasticHeavyIonKTFlux final : public KTFlux {
public:
  explicit KleinElasticHeavyIonKTFlux(const ParametersList& params)
      : KTFlux(params),
        ff_(FormFactorsFactory::get().build(formFactorsParameters(params))),
        hi_(HeavyIon::fromPdgId(ff_->pdgId())),  // the form factors (set by the beam) define the ion
        mass2_(hi_.mass() * hi_.mass()) {}

  static ParametersDescription description() {
    auto desc = KTFlux::description();
    desc.setDescription("Elastic photon emission from heavy ion (from Starlight)");
    desc.addAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon::Pb());
    desc.add("formFactors", FormFactorsFactory::get().describeParameters("HeavyIonDipole"));
    return desc;
  }

  bool fragmenting() const override { return false; }
  double mass2() const override { return mass2_; }
  spdgid_t partonPdgId() const override { return PDG::photon; }
  bool singleThreaded() const override { return KTFlux::singleThreaded() || ff_->singleThreaded(); }

  double fluxQ2(double x, double kt2, double q2) const override {
    if (!x_range_.contains(x))
      return 0.;

    const auto ff = (*ff_)(q2);
    const double ela1 = pow(kt2 / q2 / (1. - x), 2);
    const double ela2 = pow(ff.GE, 2);
    //const double ela3 = kt2 / q2;
    const auto z = static_cast<unsigned short>(hi_.Z);
    return alpha_over_pi_ * z * z * ela1 * ela2 / q2;
  }

  double fluxMX2(double x, double kt2, double) const override { return fluxQ2(x, kt2, utils::kt::q2(x, kt2, mass2())); }

private:
  /// Form factors modelling parameters, defaulting to the flux-defined heavy ion if no ion was set by the beam
  static ParametersList formFactorsParameters(const ParametersList& params) {
    auto ff_params = params.get<ParametersList>("formFactors");
    if (!ff_params.has<int>("pdgId"))
      ff_params.set<int>("pdgId", static_cast<int>(params.get<pdgid_t>("heavyIon", static_cast<pdgid_t>(HeavyIon::Pb()))));
    return ff_params;
  }
  const std::unique_ptr<formfac::Parameterisation> ff_;
  const HeavyIon hi_;
  const double mass2_;
};
REGISTER_KT_FLUX("KleinElasticHI", 101, KleinElasticHeavyIonKTFlux);
