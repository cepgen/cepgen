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
#include "CepGenAddOns/CubaWrapper/IntegratorCuba.h"

namespace cepgen {
  static Integrand* gIntegrand{nullptr};

  IntegratorCuba::IntegratorCuba(const ParametersList& params)
      : Integrator(params),
        ncomp_(steer<int>("ncomp")),
        nvec_(steer<int>("nvec")),
        epsrel_(steer<double>("epsrel")),
        epsabs_(steer<double>("epsabs")),
        mineval_(steer<int>("mineval")),
        maxeval_(steer<int>("maxeval")),
        verbose_(steer<int>("verbose")) {}

  void IntegratorCuba::setIntegrand(Integrand& integr) {
    Integrator::setIntegrand(integr);
    gIntegrand = integrand_;
  }

  ParametersDescription IntegratorCuba::description() {
    auto desc = Integrator::description();
    desc.setDescription("Cuba generic integration algorithm");
    desc.add<int>("ncomp", 1).setDescription("number of components of the integrand");
    desc.add<int>("nvec", 1).setDescription("number of samples received by the integrand");
    desc.add<double>("epsrel", 1.e-3).setDescription("requested relative accuracy");
    desc.add<double>("epsabs", 1.e-12).setDescription("requested absolute accuracy");
    desc.add<int>("mineval", 0).setDescription("minimum number of integrand evaluations required");
    desc.add<int>("maxeval", 50000).setDescription("(approximate) maximum number of integrand evaluations allowed");
    desc.add<int>("verbose", 0);
    return desc;
  }

  int cuba_integrand(const int* ndim, const double xx[], const int* /*ncomp*/, double ff[], void* /*userdata*/) {
    if (!gIntegrand)
      throw CG_FATAL("cuba_integrand") << "Integrand not set for the Cuba algorithm!";
    //FIXME handle the non-[0,1] ranges
    ff[0] = gIntegrand->eval(std::vector<double>(xx, xx + *ndim));
    return 0;
  }
}  // namespace cepgen
