/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2023  Laurent Forthomme
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

#include "CepGen/CollinearFluxes/Parameterisation.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/FunctionsWrappers.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  namespace collflux {
    class GammaIntegrated : public Parameterisation {
    public:
      explicit GammaIntegrated(const ParametersList& params)
          : Parameterisation(params),
            flux_(PartonFluxFactory::get().build(steer<ParametersList>("ktPartonFlux"))),
            integr_(AnalyticIntegratorFactory::get().build(params.get<ParametersList>("analyticalIntegrator"))) {
        if (!flux_->ktFactorised())
          throw CG_FATAL("GammaIntegrated") << "Input flux has to be unintegrated.";
        // initialise the function to integrate
        func_.reset(new utils::Function1D([&](double kt2, void* params) {
          const auto& args = *static_cast<FluxArguments*>(params);
          return (*flux_)(args.x, kt2, args.mf2) / kt2;
        }));
        CG_INFO("GammaIntegrated") << "kt flux-integrated collinear flux evaluator initialised.\n\t"
                                   << "Q^2 integration range: " << q2_range_ << " GeV^2\n\t"
                                   << "Unintegrated flux: " << flux_->description() << ".";
      }

      bool fragmenting() const override final { return flux_->fragmenting(); }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("kt-integrated photon flux");
        desc.addAs<ParametersDescription>("ktPartonFlux",
                                          //PartonFluxFactory::get().describeParameters("BudnevElasticKT"))
                                          ParametersDescription().setName<std::string>("BudnevElasticKT"))
            .setDescription("Type of unintegrated kT-dependent parton flux");
        desc.add<ParametersDescription>("analyticalIntegrator", ParametersDescription().setName<std::string>("gsl"))
            .setDescription("Steering parameters for the analytical integrator");
        CG_LOG << desc;
        return desc;
      }

      double operator()(double x, double mx2) const override {
        static const Limits x_valid_range{0., 1.};
        if (x == 0. || !x_valid_range.contains(x))
          return 0.;
        const FluxArguments params{x, mx2};
        return 2. * M_PI * integr_->integrate(*func_, params, q2_range_) / x;
      }

    private:
      std::unique_ptr<PartonFlux> flux_;
      std::unique_ptr<utils::Function1D> func_;
      std::unique_ptr<AnalyticIntegrator> integr_;

      struct FluxArguments {
        double x{0.}, mf2{0.};
      };
    };
  }  // namespace collflux
}  // namespace cepgen
typedef cepgen::collflux::GammaIntegrated CF_GI;
REGISTER_COLLFLUX("GammaIntegrated", CF_GI)
