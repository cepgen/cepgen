/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#include "CepGen/Integration/FunctionIntegrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Utils/FunctionWrapper.h"

using namespace cepgen;

Integrator::Integrator(const ParametersList& params)
    : NamedModule(params),
      integrand_parameters_(steer<ParametersList>("integrandParameters")),
      verbosity_(steer<int>("verbosity")) {}

double Integrator::eval(Integrand& integrand, const std::vector<double>& x) const { return integrand.eval(x); }

Value Integrator::integrate(Integrand& integrand, const std::vector<Limits>& range) {
  if (range.size() < integrand.size()) {
    auto normalised_range = range;
    for (size_t i = 0; i < integrand.size() - range.size(); ++i)
      normalised_range.emplace_back(0., 1.);
    return run(integrand, normalised_range);
  }
  return run(integrand, range);
}

Value Integrator::integrate(const std::function<double(double)>& integrand, const Limits& range_1d) {
  auto function_integrand =
      FunctionIntegrand{1, [&integrand](const std::vector<double>& x) { return integrand(x.at(0)); }};
  return integrate(function_integrand, std::vector{range_1d});
}

Value Integrator::integrate(const std::function<double(const std::vector<double>&)>& integrand,
                            const std::vector<Limits>& range) {
  auto function_integrand = FunctionIntegrand{range.size(), integrand};
  return integrate(function_integrand, range);
}

ParametersDescription Integrator::description() {
  auto desc = ParametersDescription();
  desc.add("integrandParameters", ParametersDescription{}).setDescription("parameters for the integrand");
  desc.add("verbosity", 0).setDescription("integrator verbosity");
  return desc;
}
