/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

#include <cuba.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGenAddOns/CubaWrapper/CubaIntegrator.h"

namespace cepgen {
  /// Cuba implementation of the VEGAS integration algorithm
  class CuhreCubaIntegrator final : public CubaIntegrator {
  public:
    explicit CuhreCubaIntegrator(const ParametersList& params) : CubaIntegrator(params), key_(steer<int>("key")) {
      CG_DEBUG("Integrator:build") << "Cuba-Cuhre integrator built.";
    }

    static ParametersDescription description() {
      auto desc = CubaIntegrator::description();
      desc.setDescription("Cuba implementation of the Cuhre algorithm");
      desc.add<int>("key", 0).setDescription(
          "basic integration rule:\n"
          "key = 7, 9, 11, 13 selects the cubature rule of degree key. Note that the degree-11\n"
          "rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.\n"
          "For other values, the default rule is taken, which is the degree-13 rule in 2 dimensions,\n"
          "the degree-11 rule in 3 dimensions, and the degree-9 rule otherwise.");
      return desc;
    }

    Value integrate() override {
      int nregions, neval, fail;
      double integral, error, prob;

      Cuhre(gIntegrand->size(),
            ncomp_,
            cuba_integrand,
            nullptr,
            nvec_,
            epsrel_,
            epsabs_,
            verbose_,
            mineval_,
            maxeval_,
            key_,
            nullptr,  // statefile
            nullptr,  // spin
            &nregions,
            &neval,
            &fail,
            &integral,
            &error,
            &prob);

      CG_DEBUG("CuhreCubaIntegrator:integrate")
          << "Number of regions needed: " << nregions << ".\nNumber of function evaluations: " << neval
          << "\nError flag: " << fail << ".";

      return Value{integral, error};
    }

  private:
    const int key_;
  };
}  // namespace cepgen
REGISTER_INTEGRATOR("cuba-cuhre", CuhreCubaIntegrator);
