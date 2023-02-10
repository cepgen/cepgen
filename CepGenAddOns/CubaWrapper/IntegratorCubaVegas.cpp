/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2022  Laurent Forthomme
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
#include "CepGenAddOns/CubaWrapper/IntegratorCuba.h"
#include "cuba.h"

namespace cepgen {
  /// Cuba implementation of the VEGAS integration algorithm
  class IntegratorCubaVegas : public IntegratorCuba {
  public:
    explicit IntegratorCubaVegas(const ParametersList&);
    void integrate(double&, double&) override;

    static ParametersDescription description();

  private:
    int nstart_, nincrease_, nbatch_;
    int gridno_;
  };

  IntegratorCubaVegas::IntegratorCubaVegas(const ParametersList& params)
      : IntegratorCuba(params),
        nstart_(steer<int>("NStart")),
        nincrease_(steer<int>("NIncrease")),
        nbatch_(steer<int>("NBatch")),
        gridno_(steer<int>("GridNo")) {
    //--- a bit of printout for debugging
    CG_DEBUG("Integrator:build") << "Cuba-VEGAS integrator built.";
  }

  void IntegratorCubaVegas::integrate(double& result, double& abserr) {
    int neval, fail;
    double integral, error, prob;

    Vegas(integrand_->size(),
          ncomp_,
          cuba_integrand,
          nullptr,
          nvec_,
          epsrel_,
          epsabs_,
          verbose_,
          seed_,
          mineval_,
          maxeval_,
          nstart_,
          nincrease_,
          nbatch_,
          gridno_,
          nullptr,
          nullptr,
          &neval,
          &fail,
          &integral,
          &error,
          &prob);

    result_ = result = integral;
    err_result_ = abserr = error;
  }

  ParametersDescription IntegratorCubaVegas::description() {
    auto desc = IntegratorCuba::description();
    desc.setDescription("Cuba implementation of the VEGAS algorithm");
    desc.add<int>("NStart", 1000).setDescription("number of integrand evaluations per iteration to start with");
    desc.add<int>("NIncrease", 500).setDescription("increase in the number of integrand evaluations per iteration");
    desc.add<int>("NBatch", 1000)
        .setDescription("number of points sent in one MathLink packet to be sampled by Mathematica");
    desc.add<int>("GridNo", 0).setDescription("slot in the internal grid table");
    return desc;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("cuba-vegas", IntegratorCubaVegas);
