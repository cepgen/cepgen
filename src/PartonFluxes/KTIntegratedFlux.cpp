/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2025  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/CollinearFlux.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/FunctionWrapper.h"
#include "CepGen/Utils/Limits.h"

using namespace cepgen;

class KTIntegratedFlux : public CollinearFlux {
public:
  explicit KTIntegratedFlux(const ParametersList& params)
      : CollinearFlux(params),
        integrator_(IntegratorFactory::get().build(steer<ParametersList>("integrator"))),
        flux_(KTFluxFactory::get().build(steer<ParametersList>("ktFlux"))),
        kt2_range_(steer<Limits>("kt2range")) {
    if (!flux_->ktFactorised())
      throw CG_FATAL("GammaIntegrated") << "Input flux has to be unintegrated.";
    // initialise the functions to integrate
    CG_INFO("KTIntegratedFlux") << "kt flux-integrated collinear flux evaluator initialised.\n\t"
                                << "Integrator: " << integrator_->name() << "\n\t"
                                << "Q^2 integration range: " << kt2_range_ << " GeV^2\n\t"
                                << "Unintegrated flux: " << flux_->name() << ".";
  }

  bool fragmenting() const final { return flux_->fragmenting(); }
  spdgid_t partonPdgId() const final { return flux_->partonPdgId(); }
  double mass2() const final { return flux_->mass2(); }

  static ParametersDescription description() {
    auto desc = CollinearFlux::description();
    desc.setDescription("kt-integrated coll.flux");
    desc.add("integrator", IntegratorFactory::get().describeParameters("gsl"))
        .setDescription("Steering parameters for the analytical integrator");
    desc.add("ktFlux", PartonFluxFactory::get().describeParameters("BudnevElastic"))
        .setDescription("Type of unintegrated kT-dependent parton flux");
    desc.add("kt2range", Limits{0., 1.e4})
        .setDescription("kinematic range for the parton transverse virtuality, in GeV^2");
    return desc;
  }

  double fluxQ2(double x, double q2) const override {
    if (!x_range_.contains(x, true))
      return 0.;
    return 2. * M_PI *
           integrator_->integrate([this, &x, &q2](double kt2) { return flux_->fluxQ2(x, kt2, q2); }, kt2_range_);
  }

  double fluxMX2(double x, double mx2) const override {
    if (!x_range_.contains(x, true))
      return 0.;
    return 2. * M_PI *
           integrator_->integrate([this, &x, &mx2](double kt2) { return flux_->fluxMX2(x, kt2, mx2); }, kt2_range_);
  }

private:
  const std::unique_ptr<Integrator> integrator_;
  const std::unique_ptr<KTFlux> flux_;
  const Limits kt2_range_;
};
REGISTER_COLLINEAR_FLUX("KTIntegrated", KTIntegratedFlux);
