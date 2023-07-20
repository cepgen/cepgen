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

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  class ElasticNucleonKTFlux : public KTFlux {
  public:
    explicit ElasticNucleonKTFlux(const ParametersList& params)
        : KTFlux(params), ff_(FormFactorsFactory::get().build(params_ + steer<ParametersList>("formFactors"))) {
      if (!ff_)
        throw CG_FATAL("ElasticNucleonKTFlux")
            << "Elastic kT flux requires a modelling of electromagnetic form factors!";
    }

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Nucl. el. photon emission");
      desc.add<ParametersDescription>("formFactors", ParametersDescription().setName<std::string>("StandardDipole"));
      return desc;
    }
    bool fragmenting() const override final { return false; }
    double mass2() const override { return mp2_; }
    pdgid_t partonPdgId() const override { return PDG::photon; }

    double operator()(double x, double kt2, double) const override {
      if (!x_range_.contains(x))
        return 0.;
      const auto q2vals = computeQ2(x, kt2);
      const double qnorm = 1. - q2vals.min / q2vals.q2;
      const auto& formfac = (*ff_)(q2vals.q2);
      return prefactor_ * formfac.FE * qnorm * qnorm / q2vals.q2;
    }

  protected:
    /// Elastic form factors computation
    std::unique_ptr<formfac::Parameterisation> ff_;
  };

  struct BudnevElasticNucleonKTFlux : public ElasticNucleonKTFlux {
    using ElasticNucleonKTFlux::ElasticNucleonKTFlux;
    static ParametersDescription description() {
      auto desc = ElasticNucleonKTFlux::description();
      desc.setDescription("Nucl. el. photon emission (Budnev flux)");
      return desc;
    }
    double operator()(double x, double kt2, double) const override final {
      if (!x_range_.contains(x))
        return 0.;
      const auto q2vals = computeQ2(x, kt2);
      const double qnorm = 1. - q2vals.min / q2vals.q2;
      const auto& formfac = (*ff_)(q2vals.q2);
      const double f_D = formfac.FE * (1. - x) * qnorm;
      const double f_C = formfac.FM;
      return prefactor_ * (f_D + 0.5 * x * x * f_C) * (1. - x) / q2vals.q2;
    }
  };

  class BudnevElasticElectronKTFlux final : public BudnevElasticNucleonKTFlux {
  public:
    explicit BudnevElasticElectronKTFlux(const ParametersList& params)
        : BudnevElasticNucleonKTFlux(params), me2_(std::pow(PDG::get().mass(PDG::electron), 2)) {}
    static ParametersDescription description() {
      auto desc = BudnevElasticNucleonKTFlux::description();
      desc.setDescription("Electron el. photon emission (Budnev flux)");
      desc.add<ParametersDescription>("formFactors", ParametersDescription().setName<std::string>("PointLikeFermion"));
      return desc;
    }
    double mass2() const override { return me2_; }

  private:
    const double me2_;
  };

  class ElasticHeavyIonKTFlux final : public ElasticNucleonKTFlux {
  public:
    explicit ElasticHeavyIonKTFlux(const ParametersList& params)
        : ElasticNucleonKTFlux(params),
          hi_(HeavyIon::fromPdgId(steer<pdgid_t>("heavyIon"))),
          mass2_(hi_.mass() * hi_.mass()) {}

    static ParametersDescription description() {
      auto desc = ElasticNucleonKTFlux::description();
      desc.setDescription("HI el. photon emission");
      desc.addAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon::Pb());
      desc.add<ParametersDescription>("formFactors", ParametersDescription().setName<std::string>("HeavyIonDipole"));
      return desc;
    }

    double mass2() const override { return mass2_; }
    double operator()(double x, double kt2, double mx2) const override {
      const auto z = (unsigned short)hi_.Z;
      return z * z * ElasticNucleonKTFlux::operator()(x, kt2, mx2);
    }

  private:
    const HeavyIon hi_;
    const double mass2_;
  };
}  // namespace cepgen

REGISTER_FLUX("ElasticKT", ElasticNucleonKTFlux);
REGISTER_FLUX("BudnevElasticKT", BudnevElasticNucleonKTFlux);
REGISTER_FLUX("BudnevElasticElectronKT", BudnevElasticElectronKTFlux);
REGISTER_FLUX("ElasticHeavyIonKT", ElasticHeavyIonKTFlux);
