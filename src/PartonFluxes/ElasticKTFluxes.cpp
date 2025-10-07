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

/// Base class for all coherent, elastic kt-dependent photon emission from beam modellings
class ElasticKTFlux : public KTFlux {
public:
  explicit ElasticKTFlux(const ParametersList& params)
      : KTFlux(params), form_factors_(FormFactorsFactory::get().build(steer<ParametersList>("formFactors"))) {
    if (!form_factors_)
      throw CG_FATAL("ElasticKTFlux") << "Elastic kT flux requires a modelling of electromagnetic form factors!";
  }

  static ParametersDescription description() {
    auto desc = KTFlux::description();
    desc.setDescription("Elastic photon emission");
    desc.add("formFactors", FormFactorsFactory::get().describeParameters("StandardDipole"));
    return desc;
  }
  bool fragmenting() const final { return false; }
  double mass2() const override { return mp2_; }
  spdgid_t partonPdgId() const override { return PDG::photon; }

  double fluxMX2(double x, double kt2, double) const override {
    if (!x_range_.contains(x))
      return 0.;
    const auto q2 = utils::kt::q2(x, kt2, mass2()), q2min = q2 - kt2 / (1. - x);
    const double q_norm = 1. - q2min / q2;
    const auto& form_fac = (*form_factors_)(q2);
    return alpha_over_pi_ * form_fac.FE * q_norm * q_norm / q2;
  }

protected:
  const std::unique_ptr<formfac::Parameterisation> form_factors_;  ///< elastic form factors modelling
};

/// Budnev coherent photon emission from a beam particle
struct BudnevElasticKTFlux : ElasticKTFlux {
  using ElasticKTFlux::ElasticKTFlux;
  static ParametersDescription description() {
    auto desc = ElasticKTFlux::description();
    desc.setDescription("Elastic photon emission (Budnev)");
    return desc;
  }
  double fluxMX2(double x, double kt2, double) const final {
    if (!x_range_.contains(x))
      return 0.;
    const auto q2 = utils::kt::q2(x, kt2, mass2()), q2min = q2 - kt2 / (1. - x), q_norm = 1. - q2min / q2;
    const auto& form_fac = (*form_factors_)(q2);
    const double f_D = form_fac.FE * (1. - x) * q_norm;
    const double f_C = form_fac.FM;
    return alpha_over_pi_ * (f_D + 0.5 * x * x * f_C) * (1. - x) / q2;
  }
};

/// Budnev coherent photon emission from a lepton beam
class BudnevElasticLeptonKTFlux final : public BudnevElasticKTFlux {
public:
  explicit BudnevElasticLeptonKTFlux(const ParametersList& params)
      : BudnevElasticKTFlux(params), ml2_(std::pow(PDG::get().mass(form_factors_->pdgId()), 2)) {
    CG_DEBUG("BudnevElasticLeptonKTFlux")
        << "Elastic kt-dependent parton-from-lepton initialised. Lepton: " << form_factors_->pdgId()
        << " (m=" << std::sqrt(ml2_) << " GeV).";
  }
  static ParametersDescription description() {
    auto desc = BudnevElasticKTFlux::description();
    desc.setDescription("Lepton elastic photon emission (Budnev)");
    desc.add("formFactors", FormFactorsFactory::get().describeParameters("PointLikeFermion"));
    return desc;
  }
  double mass2() const override { return ml2_; }

private:
  const double ml2_;
};

/// Photon emission from heavy ion
class ElasticHeavyIonKTFlux final : public ElasticKTFlux {
public:
  explicit ElasticHeavyIonKTFlux(const ParametersList& params)
      : ElasticKTFlux(params), hi_(HeavyIon::fromPdgId(form_factors_->pdgId())), mass2_(hi_.mass() * hi_.mass()) {
    CG_DEBUG("ElasticHeavyIonKTFlux") << "KT-factorised elastic photon-from-HI flux evaluator built for HI=" << hi_
                                      << ", (mass=" << hi_.mass()
                                      << "), electromagnetic form factors: " << form_factors_->parameters() << ".";
  }

  static ParametersDescription description() {
    auto desc = ElasticKTFlux::description();
    desc.setDescription("HI elastic photon emission");
    desc.add("formFactors", FormFactorsFactory::get().describeParameters("HeavyIonDipole"));
    return desc;
  }

  double mass2() const override { return mass2_; }

  double fluxMX2(double x, double kt2, double mx2) const override {
    const auto z = static_cast<unsigned short>(hi_.Z);
    return z * z * ElasticKTFlux::fluxMX2(x, kt2, mx2);
  }

private:
  const HeavyIon hi_;
  const double mass2_;
};
REGISTER_KT_FLUX("Elastic", 0, ElasticKTFlux);
REGISTER_KT_FLUX("BudnevElastic", 10, BudnevElasticKTFlux);
REGISTER_KT_FLUX("BudnevElasticLepton", 12, BudnevElasticLeptonKTFlux);
REGISTER_KT_FLUX("ElasticHeavyIon", 100, ElasticHeavyIonKTFlux);
