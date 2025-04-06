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

#include <cmath>

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Message.h"

namespace cepgen::formfac {
  class InelasticNucleon : public Parameterisation {
  public:
    explicit InelasticNucleon(const ParametersList& params)
        : Parameterisation(params),
          sf_(StructureFunctionsFactory::get().build(steer<ParametersList>("structureFunctions"))),
          integrator_(IntegratorFactory::get().build(steer<ParametersList>("integrator"))),
          compute_fm_(steer<bool>("computeFM")),
          mx_range_(steer<Limits>("mxRange")),
          mx2_range_{mx_range_.min() * mx_range_.min(), mx_range_.max() * mx_range_.max()},
          dm2_range_{mx2_range_.min() - mp2_, mx2_range_.max() - mp2_},
          eval_fe_([this](double mx2) {
            const auto xbj = utils::xBj(q2_, mp2_, mx2);
            return sf_->F2(xbj, q2_) * xbj;
          }),
          eval_fm_([this](double mx2) {
            const auto xbj = utils::xBj(q2_, mp2_, mx2);
            return sf_->F2(xbj, q2_) / xbj;
          }) {
      CG_INFO("InelasticNucleon") << "Inelastic nucleon form factors parameterisation built with:\n"
                                  << " * structure functions modelling: " << steer<ParametersList>("structureFunctions")
                                  << "\n"
                                  << " * integrator algorithm: " << steer<ParametersList>("integrator") << "\n"
                                  << " * diffractive mass range: " << steer<Limits>("mxRange") << " GeV^2.";
    }

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Proton inelastic (SF)");
      desc.add("structureFunctions", StructureFunctionsFactory::get().describeParameters("LUXLike"))
          .setDescription("type of structure functions parameterisation for the dissociative emission");
      desc.add("integrator", IntegratorFactory::get().describeParameters("gsl"))
          .setDescription("type of numerical integrator algorithm to use");
      desc.add("computeFM", false).setDescription("compute, or neglect the F2/xbj^3 term");
      desc.add("mxRange", Limits{1.0732 /* mp + mpi0 */, 20.}).setDescription("diffractive mass range (in GeV/c^2)");
      return desc;
    }

  protected:
    void eval() override {
      const auto inv_q2 = 1. / q2_;
      const auto fe = integrator_->integrate(eval_fe_, mx2_range_) * inv_q2;
      const auto fm = compute_fm_ ? integrator_->integrate(eval_fm_, mx2_range_) * inv_q2 : 0.;
      setFEFM(fe, fm);
    }
    bool fragmenting() const override { return true; }

  private:
    const std::unique_ptr<strfun::Parameterisation> sf_;
    const std::unique_ptr<Integrator> integrator_;
    const double compute_fm_;
    const Limits mx_range_, mx2_range_, dm2_range_;
    const std::function<double(double)> eval_fe_, eval_fm_;
  };
}  // namespace cepgen::formfac
using cepgen::formfac::InelasticNucleon;
REGISTER_FORMFACTORS("InelasticNucleon", InelasticNucleon);
