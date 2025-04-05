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

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"

using namespace cepgen;

/// Base class for all coherent, elastic kt-dependent photon emission from nucleons modellings
class ElasticNucleonKTFlux : public KTFlux {
public:
  explicit ElasticNucleonKTFlux(const ParametersList& params)
      : KTFlux(params), form_factors_(FormFactorsFactory::get().build(steer<ParametersList>("formFactors"))) {
    if (!form_factors_)
      throw CG_FATAL("ElasticNucleonKTFlux") << "Elastic kT flux requires a modelling of electromagnetic form factors!";
  }

  static ParametersDescription description() {
    auto desc = KTFlux::description();
    desc.setDescription("Nucl. el. photon emission");
    desc.add("formFactors", FormFactorsFactory::get().describeParameters("StandardDipole"));
    return desc;
  }
  bool fragmenting() const final { return false; }
  double mass2() const override { return mp2_; }
  pdgid_t partonPdgId() const override { return PDG::photon; }

  double fluxMX2(double x, double kt2, double) const override {
    if (!x_range_.contains(x))
      return 0.;
    const auto q2 = utils::kt::q2(x, kt2, mass2()), q2min = q2 - kt2 / (1. - x);
    const double qnorm = 1. - q2min / q2;
    const auto& formfac = (*form_factors_)(q2);
    return alpha_over_pi_ * formfac.FE * qnorm * qnorm / q2;
  }

protected:
  const std::unique_ptr<formfac::Parameterisation> form_factors_;  ///< elastic form factors modelling
};

/// Budnev coherent photon emission from a nucleon
struct BudnevElasticNucleonKTFlux : ElasticNucleonKTFlux {
  using ElasticNucleonKTFlux::ElasticNucleonKTFlux;
  static ParametersDescription description() {
    auto desc = ElasticNucleonKTFlux::description();
    desc.setDescription("Nucl. el. photon emission (Budnev flux)");
    return desc;
  }
  double fluxMX2(double x, double kt2, double) const final {
    if (!x_range_.contains(x))
      return 0.;
    const auto q2 = utils::kt::q2(x, kt2, mass2()), q2min = q2 - kt2 / (1. - x), qnorm = 1. - q2min / q2;
    const auto& ff_value = (*form_factors_)(q2);
    const double f_D = ff_value.FE * (1. - x) * qnorm;
    const double f_C = ff_value.FM;
    return alpha_over_pi_ * (f_D + 0.5 * x * x * f_C) * (1. - x) / q2;
  }
};

/// Budnev coherent photon emission from a lepton beam
class BudnevElasticLeptonKTFlux final : public BudnevElasticNucleonKTFlux {
public:
  explicit BudnevElasticLeptonKTFlux(const ParametersList& params)
      : BudnevElasticNucleonKTFlux(params), ml2_(std::pow(steer<ParticleProperties>("pdgId").mass, 2)) {}
  static ParametersDescription description() {
    auto desc = BudnevElasticNucleonKTFlux::description();
    desc.setDescription("Lepton el. photon emission (Budnev flux)");
    desc.add("formFactors", FormFactorsFactory::get().describeParameters("PointLikeFermion"));
    desc.addAs<pdgid_t>("pdgId", PDG::electron).setDescription("lepton flavour");
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
                                      << "), electromagnetic form factors: " << form_factors_->parameters() << ".";
  }

  static ParametersDescription description() {
    auto desc = ElasticNucleonKTFlux::description();
    desc.setDescription("HI el. photon emission");
    desc.addAs<pdgid_t>("heavyIon", HeavyIon::Pb());
    desc.add("formFactors", FormFactorsFactory::get().describeParameters("HeavyIonDipole"));
    return desc;
  }

  double mass2() const override { return mass2_; }

  double fluxMX2(double x, double kt2, double mx2) const override {
    const auto z = static_cast<unsigned short>(hi_.Z);
    return z * z * ElasticNucleonKTFlux::fluxMX2(x, kt2, mx2);
  }

private:
  const HeavyIon hi_;
  const double mass2_;
};
REGISTER_KT_FLUX("Elastic", 0, ElasticNucleonKTFlux);
REGISTER_KT_FLUX("BudnevElastic", 10, BudnevElasticNucleonKTFlux);
REGISTER_KT_FLUX("BudnevElasticLepton", 12, BudnevElasticLeptonKTFlux);
REGISTER_KT_FLUX("ElasticHeavyIon", 100, ElasticHeavyIonKTFlux);
