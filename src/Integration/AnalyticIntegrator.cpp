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

#include "CepGen/Integration/AnalyticIntegrator.h"

#include "CepGen/Utils/FunctionWrapper.h"

using namespace cepgen;

AnalyticIntegrator::AnalyticIntegrator(const ParametersList& params)
    : NamedModule(params),
      range_(steer<Limits>("range")),
      func_params_(steer<ParametersList>("params")),
      verbosity_(steer<int>("verbosity")) {}

double AnalyticIntegrator::integrate(const std::function<double(double)>& func, const Limits& lim) const {
  return run(utils::FunctionWrapper(func), nullptr, lim);
}

ParametersDescription AnalyticIntegrator::description() {
  auto desc = ParametersDescription();
  desc.add<Limits>("range", Limits{0., 1.}).setDescription("integration range");
  desc.add<ParametersDescription>("params", ParametersDescription())
      .setDescription("parameters for the function to be integrated");
  desc.add<int>("verbosity", 0).setDescription("integrator verbosity");
  return desc;
}