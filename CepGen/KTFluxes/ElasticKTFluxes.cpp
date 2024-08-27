/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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
#include "CepGen/Physics/Utils.h"

namespace cepgen {
  /// Base class for all coherent, elastic kt-dependent photon emission from nucleons modellings
  class ElasticNucleonKTFlux : public KTFlux {
  public:
    explicit ElasticNucleonKTFlux(const ParametersList& params)
        : KTFlux(params), ff_(FormFactorsFactory::get().build(steer<ParametersList>("formFactors"))) {
      if (!ff_)
        throw CG_FATAL("ElasticNucleonKTFlux")
            << "Elastic kT flux requires a modelling of electromagnetic form factors!";
    }

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Nucl. el. photon emission");
      desc.add<ParametersDescription>("formFactors", FormFactorsFactory::get().describeParameters("StandardDipole"));
      return desc;
    }
    bool fragmenting() const override final { return false; }
    double mass2() const override { return mp2_; }
    pdgid_t partonPdgId() const override { return PDG::photon; }

    double fluxQ2(double x, double kt2, double q2) const override {
      if (!x_range_.contains(x))
        return 0.;
      const auto q2min = q2 - kt2 / (1. - x), qnorm = 1. - q2min / q2;
      const auto& formfac = (*ff_)(q2);
      return prefactor_ * formfac.FE * qnorm * qnorm / q2;
    }

  protected:
    const std::unique_ptr<formfac::Parameterisation> ff_;  ///< elastic form factors modelling
  };

  /// Budnev coherent photon emission from a nucleon
  struct BudnevElasticNucleonKTFlux : public ElasticNucleonKTFlux {
    using ElasticNucleonKTFlux::ElasticNucleonKTFlux;
    static ParametersDescription description() {
      auto desc = ElasticNucleonKTFlux::description();
      desc.setDescription("Nucl. el. photon emission (Budnev flux)");
      return desc;
    }
    double fluxQ2(double x, double kt2, double q2) const override final {
      if (!x_range_.contains(x))
        return 0.;
      const auto q2min = q2 - kt2 / (1. - x), qnorm = 1. - q2min / q2;
      const auto& formfac = (*ff_)(q2);
      const double f_D = formfac.FE * (1. - x) * qnorm;
      const double f_C = formfac.FM;
      return prefactor_ * (f_D + 0.5 * x * x * f_C) * (1. - x) / q2;
    }
  };

  /// Budnev coherent photon emission from a lepton beam
  class BudnevElasticLeptonKTFlux final : public BudnevElasticNucleonKTFlux {
  public:
    explicit BudnevElasticLeptonKTFlux(const ParametersList& params)
        : BudnevElasticNucleonKTFlux(params), ml2_(std::pow(PDG::get().mass(steer<pdgid_t>("pdgId")), 2)) {}
    static ParametersDescription description() {
      auto desc = BudnevElasticNucleonKTFlux::description();
      desc.setDescription("Lepton el. photon emission (Budnev flux)");
      desc.add<ParametersDescription>("formFactors", FormFactorsFactory::get().describeParameters("PointLikeFermion"));
      desc.add<pdgid_t>("pdgId", PDG::electron).setDescription("lepton flavour");
      return desc;
    }
    double mass2() const override { return ml2_; }

  private:
    const double ml2_;
  };

  /// Photon emission from heavy ion
  class ElasticHeavyIonKTFlux final : public ElasticNucleonKTFlux {
  public:
    explicit ElasticHeavyIonKTFlux(const ParametersList& params)
        : ElasticNucleonKTFlux(params),
          hi_(HeavyIon::fromPdgId(steer<pdgid_t>("heavyIon"))),
          mass2_(hi_.mass() * hi_.mass()) {
      CG_DEBUG("ElasticHeavyIonKTFlux") << "KT-factorised elastic photon-from-HI flux evaluator built for HI=" << hi_
                                        << ", (mass=" << hi_.mass()
                                        << "), electromagnetic form factors: " << ff_->parameters() << ".";
    }

    static ParametersDescription description() {
      auto desc = ElasticNucleonKTFlux::description();
      desc.setDescription("HI el. photon emission");
      desc.addAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon::Pb());
      desc.add<ParametersDescription>("formFactors", FormFactorsFactory::get().describeParameters("HeavyIonDipole"));
      return desc;
    }

    double mass2() const override { return mass2_; }

    double fluxQ2(double x, double kt2, double q2) const override {
      const auto z = static_cast<unsigned short>(hi_.Z);
      return z * z * ElasticNucleonKTFlux::fluxQ2(x, kt2, q2);
    }

  private:
    const HeavyIon hi_;
    const double mass2_;
  };
}  // namespace cepgen

REGISTER_KT_FLUX("Elastic", 0, ElasticNucleonKTFlux);
REGISTER_KT_FLUX("BudnevElastic", 10, BudnevElasticNucleonKTFlux);
REGISTER_KT_FLUX("BudnevElasticLepton", 12, BudnevElasticLeptonKTFlux);
REGISTER_KT_FLUX("ElasticHeavyIon", 100, ElasticHeavyIonKTFlux);
