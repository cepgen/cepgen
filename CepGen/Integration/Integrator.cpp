/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Integration/FunctionIntegrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"

namespace cepgen {
  Integrator::Integrator(const ParametersList& params)
      : NamedModule(params),
        rnd_gen_(RandomGeneratorFactory::get().build(
            steer<ParametersList>("randomGenerator").set<unsigned long long>("seed", steer<int>("seed")))),
        verbosity_(steer<int>("verbose")) {}

  void Integrator::checkLimits(const Integrand& integrand) {
    const auto ps_size = integrand.size();
    if (ps_size == 0)
      throw CG_FATAL("Integrator:checkLimits") << "Invalid phase space dimension for integrand: " << ps_size << ".";
    else if (limits_.empty())
      setLimits(std::vector<Limits>(ps_size, Limits{0., 1.}));
    else if (limits_.size() != ps_size) {
      CG_DEBUG("Integrator:checkLimits") << "Incompatible phase space size: prepared=" << limits_.size()
                                         << ", integrand=" << ps_size << ".";
      auto lims = limits_;
      const auto booked_size = lims.size();
      if (booked_size < ps_size)
        for (size_t i = 0; i < ps_size - booked_size; ++i)
          lims.emplace_back(0., 1.);
      else
        lims.resize(ps_size);
      setLimits(lims);
    }
  }

  double Integrator::eval(Integrand& integrand, const std::vector<double>& x) const { return integrand.eval(x); }

  double Integrator::uniform(const Limits& lim) const { return rnd_gen_->uniform(lim.min(), lim.max()); }

  Value Integrator::integrate(Integrand& integrand) {
    if (limits_.size() != integrand.size())
      limits_ = std::vector<Limits>(integrand.size(), Limits{0., 1.});
    return integrate(integrand);
  }

  Value Integrator::integrate(const std::function<double(const std::vector<double>&)>& func,
                              const ParametersList& params,
                              size_t num_vars) {
    return integrate(func, params, std::vector<Limits>(num_vars, Limits{0., 1.}));
  }

  Value Integrator::integrate(const std::function<double(const std::vector<double>&)>& func,
                              const ParametersList& params,
                              const std::vector<Limits>& limits) {
    auto integr = IntegratorFactory::get().build(params);
    integr->setLimits(limits);
    auto integrand = FunctionIntegrand(limits.size(), func);
    return integr->integrate(integrand);
  }

  ParametersDescription Integrator::description() {
    auto desc = ParametersDescription();
    desc.setDescription("Unnamed integrator");
    desc.add<int>("seed", time(nullptr)).setDescription("Random number generator seed");
    desc.add<int>("verbose", 1).setDescription("Verbosity level");
    desc.add<ParametersDescription>("randomGenerator", ParametersDescription().setName<std::string>("stl"))
        .setDescription("random number generator engine");
    return desc;
  }
}  // namespace cepgen
