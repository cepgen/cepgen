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
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  namespace formfac {
    class InelasticFormFactors : public Parameterisation {
    public:
      explicit InelasticFormFactors(const ParametersList& params)
          : Parameterisation(params),
            sf_(StructureFunctionsFactory::get().build(steer<ParametersList>("structureFunctions"))),
            integr_(AnalyticIntegratorFactory::get().build(steer<ParametersList>("integrator"))),
            xbj_range_(steer<Limits>("xbjRange")),
            xbjm3_range_(std::pow(xbj_range_.max(), -3), std::min(1.e9, std::pow(xbj_range_.min(), -3))) {}

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Proton inelastic (SF)");
        desc.add<ParametersDescription>("structureFunctions", ParametersDescription().setName<int>(11));
        desc.add<ParametersDescription>("integrator", ParametersDescription().setName<std::string>("gsl"));
        desc.add<Limits>("xbjRange", Limits{1.e-9, 1.});
        return desc;
      }

    protected:
      void eval() override {
        setFEFM(integr_->integrate([this](double xbj) { return sf_->F2(xbj, q2_) / xbj; }, xbj_range_),
                integr_->integrate(
                    [this](double xbjm3) {
                      const auto xbj = 1. / std::cbrt(xbjm3);
                      return sf_->F2(xbj, q2_) * xbj / 3.;
                    },
                    xbjm3_range_));
      }

    private:
      const std::unique_ptr<strfun::Parameterisation> sf_;
      const std::unique_ptr<AnalyticIntegrator> integr_;
      const Limits xbj_range_, xbjm3_range_;
    };
  }  // namespace formfac
}  // namespace cepgen
using cepgen::formfac::InelasticFormFactors;
REGISTER_FORMFACTORS("InelasticNucleon", InelasticFormFactors);
