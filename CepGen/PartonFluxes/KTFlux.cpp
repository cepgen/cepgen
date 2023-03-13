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
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  class NucleonKTFlux : public KTFlux {
  public:
    explicit NucleonKTFlux(const ParametersList& params) : KTFlux(params) {}

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.add<pdgid_t>("pdgId", PDG::proton).setDescription("PDG identifier of the incoming nucleon");
      return desc;
    }

  protected:
    virtual double mass2() const = 0;

    struct Q2Values {
      double min{0.}, q2{0.};
    };
    /// Compute the minimum and kT-dependent Q^2
    Q2Values computeQ2(double x, double kt2, double mx2 = 0.) const {
      Q2Values out;
      const auto mi2 = mass2();
      const auto dm2 = (mx2 == 0.) ? 0. : mx2 - mi2;
      out.min = ((x * dm2) + x * x * mi2) / (1. - x);
      out.q2 = out.min + kt2 / (1. - x);
      return out;
    }
  };

  class ElasticNucleonKTFlux : public NucleonKTFlux {
  public:
    explicit ElasticNucleonKTFlux(const ParametersList& params)
        : NucleonKTFlux(params), ff_(FormFactorsFactory::get().build(params.get<ParametersList>("formFactors"))) {
      if (!ff_)
        throw CG_FATAL("ElasticNucleonKTFlux")
            << "Elastic kT flux requires a modelling of electromagnetic form factors!";
    }

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Nucleon elastic photon emission");
      desc.add<ParametersDescription>("formFactors", ParametersDescription().setName<std::string>("StandardDipole"));
      return desc;
    }
    bool fragmenting() const override final { return false; }
    double mass2() const override { return mp2_; }

    double operator()(double x, double kt2, double) const override {
      if (!x_range_.contains(x))
        return 0.;
      const auto q2vals = computeQ2(x, kt2);
      const double qnorm = 1. - q2vals.min / q2vals.q2;
      const auto& formfac = (*ff_)(q2vals.q2);
      return constants::ALPHA_EM * M_1_PI * formfac.FE * qnorm * qnorm / q2vals.q2;
    }

  protected:
    /// Elastic form factors computation
    std::unique_ptr<formfac::Parameterisation> ff_;
  };

  class ElasticHeavyIonKTFlux : public ElasticNucleonKTFlux {
  public:
    explicit ElasticHeavyIonKTFlux(const ParametersList& params)
        : ElasticNucleonKTFlux(params),
          hi_(HeavyIon::fromPdgId(steer<pdgid_t>("heavyIon"))),
          mass2_(hi_.mass() * hi_.mass()) {}

    static ParametersDescription description() {
      auto desc = ElasticNucleonKTFlux::description();
      desc.setDescription("HI elastic photon emission (from Starlight)");
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

  struct BudnevElasticNucleonKTFlux final : public ElasticNucleonKTFlux {
    using ElasticNucleonKTFlux::ElasticNucleonKTFlux;
    static ParametersDescription description() {
      auto desc = ElasticNucleonKTFlux::description();
      desc.setDescription("Nucleon elastic photon emission (Budnev flux)");
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
      return constants::ALPHA_EM * M_1_PI * (f_D + 0.5 * x * x * f_C) * (1. - x) / q2vals.q2;
    }
  };

  class InelasticNucleonKTFlux : public NucleonKTFlux {
  public:
    explicit InelasticNucleonKTFlux(const ParametersList& params)
        : NucleonKTFlux(params),
          sf_(StructureFunctionsFactory::get().build(params.get<ParametersList>("structureFunctions"))) {
      if (!sf_)
        throw CG_FATAL("InelasticNucleonKTFlux") << "Inelastic kT flux requires a modelling of structure functions!";
    }

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Nucleon inelastic photon emission");
      desc.add<ParametersDescription>("structureFunctions", ParametersDescription().setName<int>(301));
      return desc;
    }

    double mass2() const override { return mp2_; }
    double operator()(double x, double kt2, double mx2) const override {
      if (!x_range_.contains(x))
        return 0.;
      if (mx2 < 0.)
        throw CG_FATAL("InelasticNucleonKTFlux") << "Diffractive mass squared mX^2 should be specified!";
      const auto q2vals = computeQ2(x, kt2, mx2);
      const auto xbj = utils::xBj(q2vals.q2, mass2(), mx2), qnorm = 1. - q2vals.min / q2vals.q2;
      return constants::ALPHA_EM * M_1_PI * sf_->F2(xbj, q2vals.q2) * (xbj / q2vals.q2) * qnorm * qnorm * (1. - x) /
             q2vals.q2;
    }

  protected:
    std::unique_ptr<strfun::Parameterisation> sf_;
  };

  struct BudnevInelasticNucleonKTFlux final : public InelasticNucleonKTFlux {
    using InelasticNucleonKTFlux::InelasticNucleonKTFlux;
    static ParametersDescription description() {
      auto desc = InelasticNucleonKTFlux::description();
      desc.setDescription("Nucleon inelastic photon emission (Budnev flux)");
      return desc;
    }
    double operator()(double x, double kt2, double mx2) const override {
      if (!x_range_.contains(x))
        return 0.;
      if (mx2 < 0.)
        throw CG_FATAL("InelasticNucleonKTFlux") << "Diffractive mass squared mX^2 should be specified!";
      const auto q2vals = computeQ2(x, kt2, mx2);
      const auto xbj = utils::xBj(q2vals.q2, mass2(), mx2), qnorm = 1. - q2vals.min / q2vals.q2;
      const double f_D = sf_->F2(xbj, q2vals.q2) * (xbj / q2vals.q2) * (1. - x) * qnorm;
      const double f_C = sf_->F1(xbj, q2vals.q2) * 2. / q2vals.q2;
      return constants::ALPHA_EM * M_1_PI * (f_D + 0.5 * x * x * f_C) * (1. - x) / q2vals.q2;
    }
  };

  /// Realistic nuclear form-factor as used in STARLIGHT
  /// See \cite Klein:2016yzr
  /*class ElasticHeavyIonKTFlux final : public KTFlux {
  public:
    explicit ElasticHeavyIonKTFlux(const ParametersList& params)
        : KTFlux(params),
          ff_(formfac::FormFactorsFactory::get().build(params.get<ParametersList>("formFactors"))) {}

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Elastic photon emission from heavy ion (from Starlight)");
      desc.addAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon::Pb());
      desc.add<ParametersDescription>("formFactors", ParametersDescription().setName<std::string>("HeavyIonDipole"));
      return desc;
    }

    bool fragmenting() const override final { return false; }

    double operator()(double x, double kt2, double) const override final {
      if (!x_range_.contains(x))
        return 0.;
      const double r_a = 1.1 * cbrt(hi_.A), a0 = 0.7, m_a = hi_.A * mp_;
      const double q2_ela = (kt2 + x * x * m_a * m_a) / (1. - x), cons = sqrt(q2_ela) / (constants::GEVM1_TO_M * 1e15);
      const double tau = cons * r_a, tau1 = cons * a0;
      const double ff1 = 3. * (sin(tau) - tau * cos(tau)) / pow(tau + 1.e-10, 3);
      const double ff2 = 1. / (1. + tau1 * tau1);
      const double ela1 = pow(kt2 / (kt2 + x * x * m_a * m_a), 2);
      const double ela2 = pow(ff1 * ff2, 2);// ela3 = 1.-( q2_ela-kt2 )/q2_ela;
      const auto z = (unsigned short)hi_.Z;
      //return constants::ALPHA_EM * M_1_PI * z * z * ela1 * ela2 / q2_ela;
      return constants::ALPHA_EM * M_1_PI * z * z * ela1 / q2_ela;
    }

  private:
    std::unique_ptr<formfac::Parameterisation> ff_;
  };*/
}  // namespace cepgen

REGISTER_FLUX("ElasticKT", ElasticNucleonKTFlux);
REGISTER_FLUX("BudnevElasticKT", BudnevElasticNucleonKTFlux);
REGISTER_FLUX("ElasticHeavyIonKT", ElasticHeavyIonKTFlux);
REGISTER_FLUX("InelasticKT", InelasticNucleonKTFlux);
REGISTER_FLUX("BudnevInelasticKT", BudnevInelasticNucleonKTFlux);
