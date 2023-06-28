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

#include <cmath>

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  /// Realistic nuclear form-factor as used in STARLIGHT
  /// See \cite Klein:2016yzr
  class KleinElasticHeavyIonKTFlux final : public KTFlux {
  public:
    explicit KleinElasticHeavyIonKTFlux(const ParametersList& params)
        : KTFlux(params),
          hi_(HeavyIon::fromPdgId(steer<pdgid_t>("heavyIon"))),
          ff_(FormFactorsFactory::get().build(params.get<ParametersList>("formFactors"))) {}

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Elastic photon emission from heavy ion (from Starlight)");
      desc.addAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon::Pb());
      desc.add<ParametersDescription>("formFactors", ParametersDescription().setName<std::string>("HeavyIonDipole"));
      return desc;
    }

    bool fragmenting() const override final { return false; }
    double mass2() const override final { return hi_.A * hi_.A * mp2_; }
    int partonPdgId() const override { return PDG::photon; }

    double operator()(double x, double kt2, double) const override final {
      if (!x_range_.contains(x))
        return 0.;

      const auto q2_ela = computeQ2(x, kt2).q2;
      const auto ff = (*ff_)(q2_ela);

      const double ela1 = pow(kt2 / q2_ela / (1. - x), 2);
      const double ela2 = pow(ff.GE, 2);
      //const double ela3 = kt2 / q2_ela;
      const auto z = (unsigned short)hi_.Z;
      return constants::ALPHA_EM * M_1_PI * z * z * ela1 * ela2 / q2_ela;
    }

  private:
    const HeavyIon hi_;
    std::unique_ptr<formfac::Parameterisation> ff_;
  };
}  // namespace cepgen

REGISTER_FLUX("KleinElasticHeavyIonKT", KleinElasticHeavyIonKTFlux);
