/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2025  Laurent Forthomme
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
#include "CepGenCuba/Integrator.h"

namespace cepgen::cuba {
  Integrand* Integrator::gIntegrand = nullptr;

  Integrator::Integrator(const ParametersList& params)
      : cepgen::Integrator(params),
        ncomp_(steer<int>("ncomp")),
        nvec_(steer<int>("nvec")),
        epsrel_(steer<double>("epsrel")),
        epsabs_(steer<double>("epsabs")),
        mineval_(steer<int>("mineval")),
        maxeval_(steer<int>("maxeval")) {}

  Value Integrator::run(Integrand& integrand, const std::vector<Limits>& /*range*/) {
    gIntegrand = &integrand;
    return integrate();
  }

  ParametersDescription Integrator::description() {
    auto desc = cepgen::Integrator::description();
    desc.setDescription("Cuba generic integration algorithm");
    desc.add("ncomp", 1).setDescription("number of components of the integrand");
    desc.add("nvec", 1).setDescription("number of samples received by the integrand");
    desc.add("epsrel", 1.e-3).setDescription("requested relative accuracy");
    desc.add("epsabs", 1.e-12).setDescription("requested absolute accuracy");
    desc.add("mineval", 0).setDescription("minimum number of integrand evaluations required");
    desc.add("maxeval", 50'000).setDescription("(approximate) maximum number of integrand evaluations allowed");
    return desc;
  }

  int cuba_integrand(const int* ndim, const double xx[], const int* /*ncomp*/, double ff[], void* /*userdata*/) {
    if (!Integrator::gIntegrand)
      throw CG_FATAL("cuba_integrand") << "Integrand not set for the Cuba algorithm!";
    //TODO: handle the non-[0,1] ranges
    ff[0] = Integrator::gIntegrand->eval(std::vector<double>(xx, xx + *ndim));
    return 0;
  }
}  // namespace cepgen::cuba
