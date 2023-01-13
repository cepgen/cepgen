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
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/PartonFlux.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  KTFlux::KTFlux(const ParametersList& params) : PartonFlux(params) {}

  ParametersDescription KTFlux::description() {
    auto desc = PartonFlux::description();
    desc.setDescription("kT-factorised flux");
    return desc;
  }

  class NucleonKTFlux : public KTFlux {
  public:
    explicit NucleonKTFlux(const ParametersList& params)
        : KTFlux(params), mi_(PDG::get().mass(steer<pdgid_t>("pdgId"))), mi2_(mi_ * mi_) {}

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.add<pdgid_t>("pdgId", 2212).setDescription("PDG identifier of the incoming nucleon");
      return desc;
    }

  protected:
    struct Q2Values {
      double min, q2;
    };
    /// Compute the minimum and kT-dependent Q^2
    Q2Values computeQ2(double x, double kt2, double mx2 = 0.) const {
      Q2Values out;
      const auto dm2 = (mx2 == 0.) ? 0. : mx2 - mi2_;
      out.min = ((x * dm2) + x * x * mi2_) / (1. - x);
      out.q2 = out.min + kt2 / (1. - x);
      return out;
    }
    const double mi_, mi2_;
  };

  class ElasticNucleonKTFlux : public NucleonKTFlux {
  public:
    explicit ElasticNucleonKTFlux(const ParametersList& params)
        : NucleonKTFlux(params),
          ff_(formfac::FormFactorsFactory::get().build(params.get<ParametersList>("formFactors"))) {
      if (!ff_)
        throw CG_FATAL("ElasticNucleonKTFlux")
            << "Elastic kT flux requires a modelling of electromagnetic form factors!";
    }

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Elastic photon emission from a nucleon");
      desc.add<std::string>("formFactors", formfac::gFFStandardDipoleHandler)
          .setDescription("Beam form factors modelling");
      return desc;
    }
    bool fragmenting() const override { return false; }

    double operator()(double x, double kt2, double) override {
      const auto q2vals = computeQ2(x, kt2);
      const double qnorm = 1. - q2vals.min / q2vals.q2;
      const auto& formfac = (*ff_)(Beam::Mode::ProtonElastic, q2vals.q2);
      return constants::ALPHA_EM * M_1_PI * formfac.FE * qnorm * qnorm / q2vals.q2;
    }

  protected:
    /// Elastic form factors computation
    std::unique_ptr<formfac::Parameterisation> ff_;
  };

  struct BudnevElasticNucleonKTFlux final : public ElasticNucleonKTFlux {
    using ElasticNucleonKTFlux::ElasticNucleonKTFlux;
    static ParametersDescription description() {
      auto desc = ElasticNucleonKTFlux::description();
      desc.setDescription("Elastic photon emission from a nucleon (Budnev flux approximation)");
      return desc;
    }
    double operator()(double x, double kt2, double) override {
      const auto q2vals = computeQ2(x, kt2);
      const double qnorm = 1. - q2vals.min / q2vals.q2;
      const auto& formfac = (*ff_)(Beam::Mode::ProtonElastic, q2vals.q2);
      const double f_D = formfac.FE * (1. - x) * qnorm;
      const double f_C = formfac.FM;
      return constants::ALPHA_EM * M_1_PI * (f_D + 0.5 * x * x * f_C) * (1. - x) / q2vals.q2;
    }
  };

  class InelasticNucleonKTFlux : public NucleonKTFlux {
  public:
    explicit InelasticNucleonKTFlux(const ParametersList& params)
        : NucleonKTFlux(params),
          sf_(strfun::StructureFunctionsFactory::get().build(params.get<ParametersList>("structureFunctions"))) {
      if (!sf_)
        throw CG_FATAL("InelasticNucleonKTFlux") << "Inelastic kT flux requires a modelling of structure functions!";
    }

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Inelastic photon emission from a nucleon");
      desc.add<ParametersDescription>(
              "structureFunctions", strfun::Parameterisation::description().setName<int>(11)  // default is SY
              )
          .setDescription("Beam inelastic structure functions modelling");
      return desc;
    }

    double operator()(double x, double kt2, double mx2) override {
      if (mx2 < 0.)
        throw CG_FATAL("InelasticNucleonKTFlux") << "Diffractive mass squared mX^2 should be specified!";
      const auto q2vals = computeQ2(x, kt2, mx2);
      const double qnorm = 1. - q2vals.min / q2vals.q2;
      const double denom = 1. / (q2vals.q2 + mx2 - mi2_), xbj = denom * q2vals.q2;
      return constants::ALPHA_EM * M_1_PI * sf_->F2(xbj, q2vals.q2) * denom * qnorm * qnorm * (1. - x) / q2vals.q2;
    }

  protected:
    std::unique_ptr<strfun::Parameterisation> sf_;
  };

  struct BudnevInelasticNucleonKTFlux final : public InelasticNucleonKTFlux {
    using InelasticNucleonKTFlux::InelasticNucleonKTFlux;
    static ParametersDescription description() {
      auto desc = InelasticNucleonKTFlux::description();
      desc.setDescription("Inelastic photon emission from a nucleon (Budnev flux approximation)");
      return desc;
    }
    double operator()(double x, double kt2, double mx2) override {
      if (mx2 < 0.)
        throw CG_FATAL("InelasticNucleonKTFlux") << "Diffractive mass squared mX^2 should be specified!";
      const auto q2vals = computeQ2(x, kt2, mx2);
      const double qnorm = 1. - q2vals.min / q2vals.q2;
      const double denom = 1. / (q2vals.q2 + mx2 - mi2_), xbj = denom * q2vals.q2;
      const double f_D = sf_->F2(xbj, q2vals.q2) * denom * (1. - x) * qnorm;
      const double f_C = sf_->F1(xbj, q2vals.q2) * 2. / q2vals.q2;
      return constants::ALPHA_EM * M_1_PI * (f_D + 0.5 * x * x * f_C) * (1. - x) / q2vals.q2;
    }
  };

  struct KMRGluonKTFlux final : public KTFlux {
    using KTFlux::KTFlux;
    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Inelastic gluon emission from proton (KMR flux modelling)");
      return desc;
    }
    double operator()(double x, double kt2, double mx2) override { return kmr::GluonGrid::get()(x, kt2, mx2); }
    int partonPdgId() const override final { return 21; }
    bool fragmenting() const override { return false; }
  };

  /// Realistic nuclear form-factor as used in STARLIGHT
  /// See \cite Klein:2016yzr
  class ElasticHeavyIonKTFlux final : public KTFlux {
  public:
    explicit ElasticHeavyIonKTFlux(const ParametersList& params)
        : KTFlux(params), hi_(steerAs<pdgid_t, HeavyIon>("heavyIon")) {}

    static ParametersDescription description() {
      auto desc = KTFlux::description();
      desc.setDescription("Elastic photon emission from heavy ion (from Starlight)");
      return desc;
    }

    bool fragmenting() const override { return false; }

    double operator()(double x, double kt2, double) override {
      const double r_a = 1.1 * cbrt(hi_.A), a0 = 0.7, m_a = hi_.A * mp_;
      const double q2_ela = (kt2 + x * x * m_a * m_a) / (1. - x), cons = sqrt(q2_ela) / (constants::GEVM1_TO_M * 1e15);
      const double tau = cons * r_a, tau1 = cons * a0;
      const double ff1 = 3. * (sin(tau) - tau * cos(tau)) / pow(tau + 1.e-10, 3);
      const double ff2 = 1. / (1. + tau1 * tau1);
      const double ela1 = pow(kt2 / (kt2 + x * x * m_a * m_a), 2);
      const double ela2 = pow(ff1 * ff2, 2) /*, ela3 = 1.-( q2_ela-kt2 )/q2_ela*/;
      const unsigned int z = (unsigned short)hi_.Z;
      return constants::ALPHA_EM * M_1_PI * z * z * ela1 * ela2 / q2_ela;
    }

  private:
    const HeavyIon hi_;
  };
}  // namespace cepgen

REGISTER_FLUX("ElasticKT", ElasticNucleonKTFlux)
REGISTER_FLUX("BudnevElasticKT", BudnevElasticNucleonKTFlux)
REGISTER_FLUX("ElasticHeavyIonKT", ElasticHeavyIonKTFlux)
REGISTER_FLUX("InelasticKT", InelasticNucleonKTFlux)
REGISTER_FLUX("BudnevInelasticKT", BudnevInelasticNucleonKTFlux)
REGISTER_FLUX("KMR", KMRGluonKTFlux)
