/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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

#include "CepGen/CollinearFluxes/CollinearFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/FunctionWrapper.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  class KTIntegratedFlux : public CollinearFlux {
  public:
    explicit KTIntegratedFlux(const ParametersList& params)
        : CollinearFlux(params),
          integr_(AnalyticIntegratorFactory::get().build(steer<ParametersList>("integrator"))),
          flux_(KTFluxFactory::get().build(steer<ParametersList>("ktFlux"))),
          kt2_range_(steer<Limits>("kt2range")),
          func_q2_([this](double kt2, void* pars) {
            const auto& args = *static_cast<std::pair<double, double>*>(pars);
            return flux_->fluxQ2(args.first, kt2, args.second);
          }),
          func_mx2_([this](double kt2, void* pars) {
            const auto& args = *static_cast<std::pair<double, double>*>(pars);
            return flux_->fluxMX2(args.first, kt2, args.second);
          }) {
      if (!flux_->ktFactorised())
        throw CG_FATAL("GammaIntegrated") << "Input flux has to be unintegrated.";
      // initialise the functions to integrate
      CG_INFO("KTIntegratedFlux") << "kt flux-integrated collinear flux evaluator initialised.\n\t"
                                  << "Analytical integrator: " << integr_->name() << "\n\t"
                                  << "Q^2 integration range: " << kt2_range_ << " GeV^2\n\t"
                                  << "Unintegrated flux: " << flux_->name() << ".";
    }

    bool fragmenting() const override final { return flux_->fragmenting(); }
    pdgid_t partonPdgId() const override final { return flux_->partonPdgId(); }
    double mass2() const override final { return flux_->mass2(); }

    static ParametersDescription description() {
      auto desc = CollinearFlux::description();
      desc.setDescription("kt-integr. coll.flux");
      desc.add<ParametersDescription>("integrator", AnalyticIntegratorFactory::get().describeParameters("gsl"))
          .setDescription("Steering parameters for the analytical integrator");
      desc.add<ParametersDescription>("ktFlux", PartonFluxFactory::get().describeParameters("BudnevElastic"))
          .setDescription("Type of unintegrated kT-dependent parton flux");
      desc.add<Limits>("kt2range", {0., 1.e4})
          .setDescription("kinematic range for the parton transverse virtuality, in GeV^2");
      return desc;
    }

    double fluxQ2(double x, double q2) const override {
      if (!x_range_.contains(x, true))
        return 0.;
      return 2. * M_PI * integr_->integrate(func_q2_, std::make_pair(x, q2), kt2_range_);
    }

    double fluxMX2(double x, double mx2) const override {
      if (!x_range_.contains(x, true))
        return 0.;
      return 2. * M_PI * integr_->integrate(func_mx2_, std::make_pair(x, mx2), kt2_range_);
    }

  private:
    const std::unique_ptr<AnalyticIntegrator> integr_;
    const std::unique_ptr<KTFlux> flux_;
    const Limits kt2_range_;
    const utils::FunctionWrapper func_q2_, func_mx2_;
  };
}  // namespace cepgen
REGISTER_COLLINEAR_FLUX("KTIntegrated", KTIntegratedFlux);
