/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGenAddOns/CubaWrapper/Integrator.h"
#include "cuba.h"

namespace cepgen::cuba {
  /// Cuba implementation of the Suave integration algorithm
  class SuaveIntegrator : public Integrator {
  public:
    explicit SuaveIntegrator(const ParametersList& params)
        : Integrator(params),
          nnew_(steer<int>("NNew")),
          nmin_(steer<int>("NMin")),
          flatness_(steer<double>("Flatness")) {
      CG_DEBUG("Integrator:build") << "Cuba-Suave integrator built.";
    }

    static ParametersDescription description() {
      auto desc = Integrator::description();
      desc.setDescription("Cuba implementation of the Suave algorithm");
      desc.add<int>("NNew", 1000).setDescription("number of new integrand evaluations in each subdivision");
      desc.add<int>("NMin", 2).setDescription(
          "minimum number of samples a former pass must contribute to a subregion to be considered in that region’s "
          "compound integral value");
      desc.add<double>("Flatness", 50.).setDescription("type of norm used to compute the fluctuation of a sample");
      return desc;
    }

    Value integrate() override {
      int neval, fail, nregions;
      double integral, error, prob;

      Suave(gIntegrand->size(),
            ncomp_,
            cuba_integrand,
            nullptr,
            nvec_,
            epsrel_,
            epsabs_,
            verbose_,
            rnd_gen_->parameters().get<unsigned long long>("seed"),
            mineval_,
            maxeval_,
            nnew_,
            nmin_,
            flatness_,
            nullptr,  // const char* statefile
            nullptr,  // void* spin
            &nregions,
            &neval,
            &fail,
            &integral,
            &error,
            &prob);
      return Value{integral, error};
    }

  private:
    const int nnew_, nmin_;
    const double flatness_;
  };
}  // namespace cepgen::cuba
using SuaveIntegratorCuba = cepgen::cuba::SuaveIntegrator;
REGISTER_INTEGRATOR("cuba_suave", SuaveIntegratorCuba);
