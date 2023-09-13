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

#include "CepGen/Core/Exception.h"
#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  class InelasticNucleonKTFlux : public KTFlux {
  public:
    explicit InelasticNucleonKTFlux(const ParametersList& params)
        : KTFlux(params),
          sf_(StructureFunctionsFactory::get().build(params.get<ParametersList>("structureFunctions"))) {
      if (!sf_)
        throw CG_FATAL("InelasticNucleonKTFlux") << "Inelastic kT flux requires a modelling of structure functions!";
    }

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Nucl. inel. photon emission");
      desc.add<ParametersDescription>("structureFunctions", ParametersDescription().setName<int>(301));
      return desc;
    }

    double mass2() const override { return mp2_; }
    pdgid_t partonPdgId() const override { return PDG::photon; }
    double fluxMX2(double x, double kt2, double mx2) const override {
      if (!x_range_.contains(x, true))
        return 0.;
      if (mx2 < 0.)
        throw CG_FATAL("InelasticNucleonKTFlux") << "Diffractive mass squared mX^2 should be specified!";
      const auto q2 = utils::kt::q2(x, kt2, mass2(), mx2), q2min = q2 - kt2 / (1. - x);
      const auto xbj = utils::xBj(q2, mass2(), mx2), qnorm = 1. - q2min / q2;
      return prefactor_ * sf_->F2(xbj, q2) * (xbj / q2) * qnorm * qnorm * (1. - x) / q2;
    }

  protected:
    std::unique_ptr<strfun::Parameterisation> sf_;
  };

  struct BudnevInelasticNucleonKTFlux final : public InelasticNucleonKTFlux {
    using InelasticNucleonKTFlux::InelasticNucleonKTFlux;
    static ParametersDescription description() {
      auto desc = InelasticNucleonKTFlux::description();
      desc.setDescription("Nucl. inel. photon emission (Budnev flux)");
      return desc;
    }
    double fluxMX2(double x, double kt2, double mx2) const override {
      if (!x_range_.contains(x, true))
        return 0.;
      if (mx2 < 0.)
        throw CG_FATAL("InelasticNucleonKTFlux") << "Diffractive mass squared mX^2 should be specified!";
      const auto q2 = utils::kt::q2(x, kt2, mass2(), mx2), q2min = q2 - kt2 / (1. - x);
      const auto xbj = utils::xBj(q2, mass2(), mx2), qnorm = 1. - q2min / q2;
      const double f_D = sf_->F2(xbj, q2) * (xbj / q2) * (1. - x) * qnorm;
      const double f_C = sf_->F1(xbj, q2) * 2. / q2;
      return prefactor_ * (f_D + 0.5 * x * x * f_C) * (1. - x) / q2;
    }
  };
}  // namespace cepgen

REGISTER_FLUX("InelasticKT", InelasticNucleonKTFlux);
REGISTER_FLUX("BudnevInelasticKT", BudnevInelasticNucleonKTFlux);
