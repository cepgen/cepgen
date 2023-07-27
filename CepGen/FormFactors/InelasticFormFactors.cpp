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
    class InelasticFormFactors : public Parameterisation {
    public:
      explicit InelasticFormFactors(const ParametersList& params)
          : Parameterisation(params),
            sf_(StructureFunctionsFactory::get().build(steer<ParametersList>("structureFunctions"))),
            integr_(AnalyticIntegratorFactory::get().build(steer<ParametersList>("integrator"))),
            mx_range_(steer<Limits>("mxRange")) {
        CG_INFO("InelasticFormFactors") << "Inelastic nucleon form factors parameterisation built with:\n"
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
        desc.add<Limits>("mxRange", Limits{1., 1000.}).setDescription("diffractive mass range (in GeV/c^2)");
        return desc;
      }

    protected:
      void eval() override {
        const auto xbj_range = Limits(utils::xBj(q2_, mp2_, mx_range_.max()), utils::xBj(q2_, mp2_, mx_range_.min())),
                   xbjm3_range = Limits(std::pow(xbj_range.max(), -3), std::pow(xbj_range.min(), -3));
        setFEFM(integr_->integrate([this](double xbj) { return sf_->F2(xbj, q2_) / xbj; }, xbj_range),
                integr_->integrate(
                    [this](double xbjm3) {
                      const auto xbj = 1. / std::cbrt(xbjm3);
                      return sf_->F2(xbj, q2_) * xbj / 3.;
                    },
                    xbjm3_range));
      }

    private:
      const std::unique_ptr<strfun::Parameterisation> sf_;
      const std::unique_ptr<AnalyticIntegrator> integr_;
      const Limits mx_range_;
    };
  }  // namespace formfac
}  // namespace cepgen
using cepgen::formfac::InelasticFormFactors;
REGISTER_FORMFACTORS("InelasticNucleon", InelasticFormFactors);
