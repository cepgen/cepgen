/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  namespace formfac {
    class InelasticNucleon : public Parameterisation {
    public:
      explicit InelasticNucleon(const ParametersList& params)
          : Parameterisation(params),
            sf_(StructureFunctionsFactory::get().build(steer<ParametersList>("structureFunctions"))),
            integr_(AnalyticIntegratorFactory::get().build(steer<ParametersList>("integrator"))),
            compute_fm_(steer<bool>("computeFM")),
            mx_range_(steer<Limits>("mxRange")),
            dm2_range_{std::pow(mx_range_.min(), 2) - mp2_, std::pow(mx_range_.max(), 2) - mp2_},
            eval_fe_([this](double xbj) { return sf_->F2(xbj, q2_) / xbj; }),
            eval_fm_([this](double xbjm3) {
              const auto xbj = inv_cbrt(xbjm3);
              return sf_->F2(xbj, q2_) * xbj;
            }) {
        CG_INFO("InelasticNucleon") << "Inelastic nucleon form factors parameterisation built with:\n"
                                    << " * structure functions modelling: "
                                    << steer<ParametersList>("structureFunctions") << "\n"
                                    << " * integrator algorithm: " << steer<ParametersList>("integrator") << "\n"
                                    << " * diffractive mass range: " << steer<Limits>("mxRange") << " GeV^2.";
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Proton inelastic (SF)");
        desc.add<ParametersDescription>("structureFunctions", ParametersDescription().setName<int>(301))
            .setDescription("type of structure functions parameterisation for the dissociative emission");
        desc.add<ParametersDescription>("integrator", ParametersDescription().setName<std::string>("gsl"))
            .setDescription("type of numerical integrator algorithm to use");
        desc.add<bool>("computeFM", true).setDescription("compute, or neglect the F2/xbj^3 term");
        desc.add<Limits>("mxRange", Limits{1.0732 /* mp + mpi0 */, 20.})
            .setDescription("diffractive mass range (in GeV/c^2)");
        return desc;
      }

    protected:
      void eval() override {
        const auto xbj_range = Limits{q2_ / (q2_ + dm2_range_.max()), q2_ / (q2_ + dm2_range_.min())};
        const auto fe = integr_->integrate(eval_fe_, xbj_range);
        double fm = 0.;
        if (compute_fm_) {
          const auto xbjm3_range = Limits{1. / xbj_range.max() / xbj_range.max() / xbj_range.max(),
                                          1. / xbj_range.min() / xbj_range.min() / xbj_range.min()};
          fm = kOneThird * integr_->integrate(eval_fm_, xbjm3_range);
        }
        setFEFM(fe, fm);
      }

      static constexpr float kOneThird = 0.33333333, kFourThirds = 1.3333333;
      /// Fast inverse cubic root
      static float inv_cbrt(float x) {
        float thirdx = x * kOneThird;
        union {  //get bits from ﬂoating-point number
          int ix;
          float fx;
        } z;
        z.fx = x;
        z.ix = 0x54a21d2a - z.ix / 3;                   // initial guess for inverse cube root
        float y = z.fx;                                 // convert integer type back to floating-point type
        y *= (kFourThirds - thirdx * y * y * y);        // 1st Newton’s iteration
        return y * (kFourThirds - thirdx * y * y * y);  // 2nd Newton’s iteration
      }

    private:
      const std::unique_ptr<strfun::Parameterisation> sf_;
      const std::unique_ptr<AnalyticIntegrator> integr_;
      const double compute_fm_;
      const Limits mx_range_, mx2_range_, dm2_range_;
      const std::function<double(double)> eval_fe_, eval_fm_;
    };
  }  // namespace formfac
}  // namespace cepgen
using cepgen::formfac::InelasticNucleon;
REGISTER_FORMFACTORS("InelasticNucleon", InelasticNucleon);
