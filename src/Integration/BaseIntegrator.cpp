/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#include "CepGen/Integration/BaseIntegrator.h"
#include "CepGen/Integration/FunctionIntegrand.h"
#include "CepGen/Utils/FunctionWrapper.h"

using namespace cepgen;

BaseIntegrator::BaseIntegrator(const ParametersList& params)
    : NamedModule(params),
      integrand_parameters_(steer<ParametersList>("params")),
      verbosity_(steer<int>("verbosity")) {}

double BaseIntegrator::eval(Integrand& integrand, const std::vector<double>& x) const { return integrand.eval(x); }

Value BaseIntegrator::integrate(const std::function<double(double)>& integrand, const Limits& range_1d) const {
  auto function_integrand =
      FunctionIntegrand{1, [&integrand](const std::vector<double>& x) { return integrand(x.at(0)); }};
  return run(function_integrand, std::vector{range_1d});
}

Value BaseIntegrator::integrate(const std::function<double(const std::vector<double>&)>& integrand,
                                const std::vector<Limits>& range) const {
  auto function_integrand = FunctionIntegrand{range.size(), integrand};
  return run(function_integrand, range);
}

ParametersDescription BaseIntegrator::description() {
  auto desc = ParametersDescription();
  desc.add("params", ParametersDescription{}).setDescription("parameters for the integrand");
  desc.add("verbosity", 0).setDescription("integrator verbosity");
  return desc;
}
