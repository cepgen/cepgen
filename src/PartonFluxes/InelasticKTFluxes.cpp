/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Math.h"

using namespace cepgen;

class InelasticNucleonKTFlux : public KTFlux {
public:
  explicit InelasticNucleonKTFlux(const ParametersList& params)
      : KTFlux(params), sf_(StructureFunctionsFactory::get().build(steer<ParametersList>("structureFunctions"))) {
    if (!sf_)
      throw CG_FATAL("InelasticNucleonKTFlux") << "Inelastic kT flux requires a modelling of structure functions!";
    CG_DEBUG("InelasticNucleonKTFlux") << "Inelastic KT-dependent flux initialised with '"
                                       << steer<ParametersList>("structureFunctions")
                                       << "' structure functions modelling.";
  }

  static ParametersDescription description() {
    auto desc = KTFlux::description();
    desc.setDescription("Nucl. inel. photon emission");
    desc.add("structureFunctions", StructureFunctionsFactory::get().describeParameters("LUXLike"));
    return desc;
  }

  double mass2() const override { return mp2_; }
  bool fragmenting() const final { return true; }
  pdgid_t partonPdgId() const override { return PDG::photon; }
  double fluxMX2(double x, double kt2, double mx2) const override {
    if (!x_range_.contains(x, true))
      return 0.;
    if (!utils::positive(mx2)) {
      CG_WARNING("InelasticNucleonKTFlux") << "Invalid diffractive mass squared mX^2 specified: " << mx2 << ".";
      return 0.;
    }
    const auto q2 = utils::kt::q2(x, kt2, mass2(), mx2), q2min = q2 - kt2 / (1. - x);
    const auto xbj = utils::xBj(q2, mass2(), mx2), qnorm = 1. - q2min / q2;
    return alpha_over_pi_ * sf_->F2(xbj, q2) * (xbj / q2) * qnorm * qnorm * (1. - x) / q2;
  }

protected:
  const std::unique_ptr<strfun::Parameterisation> sf_;
};

struct BudnevInelasticNucleonKTFlux final : public InelasticNucleonKTFlux {
  using InelasticNucleonKTFlux::InelasticNucleonKTFlux;
  static ParametersDescription description() {
    auto desc = InelasticNucleonKTFlux::description();
    desc.setDescription("Nucl. inel. photon emission (Budnev flux)");
    return desc;
  }
  double fluxMX2(double x, double kt2, double mx2) const override {
    if (!x_range_.contains(x, true))
      return 0.;
    if (!utils::positive(mx2)) {
      CG_WARNING("InelasticNucleonKTFlux") << "Invalid diffractive mass squared mX^2 specified: " << mx2 << ".";
      return 0.;
    }
    const auto q2 = utils::kt::q2(x, kt2, mass2(), mx2), q2min = q2 - kt2 / (1. - x);
    const auto xbj = utils::xBj(q2, mass2(), mx2), qnorm = 1. - q2min / q2;
    const double f_D = sf_->F2(xbj, q2) * (xbj / q2) * (1. - x) * qnorm;
    const double f_C = sf_->F1(xbj, q2) * 2. / q2;
    return alpha_over_pi_ * (f_D + 0.5 * x * x * f_C) * (1. - x) / q2;
  }
};
REGISTER_KT_FLUX("Inelastic", 1, InelasticNucleonKTFlux);
REGISTER_KT_FLUX("BudnevInelastic", 11, BudnevInelasticNucleonKTFlux);
