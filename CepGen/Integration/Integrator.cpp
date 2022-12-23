/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"

namespace cepgen {
  Integrator::Integrator(const ParametersList& params)
      : NamedModule(params),
        seed_(params_.get<int>("seed", time(nullptr))),
        verbosity_(steer<int>("verbose")),
        rnd_(0., 1.) {}

  void Integrator::setIntegrand(Integrand& integr) {
    integrand_ = &integr;
    if (limits_.size() != integrand_->size())
      limits_ = std::vector<Limits>(integrand_->size(), Limits{0., 1.});
    initialised_ = false;  // force the reinitialisation
  }

  //------------------------------------------------------------------------------------------------
  // helper / alias methods
  //------------------------------------------------------------------------------------------------

  size_t Integrator::size() const {
    if (!integrand_)
      throw CG_FATAL("Integrator:size") << "Trying to retrieve phase space size on an unitialised integrand!";
    return integrand_->size();
  }

  double Integrator::eval(const std::vector<double>& x) const {
    if (!integrand_)
      throw CG_FATAL("Integrator:eval") << "Trying to evaluate the weight on a phase space point "
                                        << "on an unitialised integrand!";
    return integrand_->eval(x);
  }

  double Integrator::uniform(double min, double max) const { return min + (max - min) * rnd_(rnd_gen_); }

  double Integrator::integrate() {
    double result, tmp;
    integrate(result, tmp);
    return result;
  }

  double Integrator::integrate(const std::function<double(const std::vector<double>&)>& func,
                               const ParametersList& params,
                               size_t num_vars) {
    return integrate(func, params, std::vector<Limits>(num_vars, Limits{0., 1.}));
  }

  double Integrator::integrate(const std::function<double(const std::vector<double>&)>& func,
                               const ParametersList& params,
                               const std::vector<Limits>& limits) {
    auto integr = IntegratorFactory::get().build(params);
    integr->setLimits(limits);
    auto integrand = FunctionIntegrand(limits.size(), func);
    integr->setIntegrand(integrand);
    return integr->integrate();
  }

  ParametersDescription Integrator::description() {
    auto desc = ParametersDescription();
    desc.setDescription("Unnamed integrator");
    desc.add<int>("seed", time(nullptr)).setDescription("Random number generator seed");
    desc.add<int>("verbose", 1).setDescription("Verbosity level");
    return desc;
  }
}  // namespace cepgen
