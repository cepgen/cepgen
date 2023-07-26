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
            xbj_range_(steer<Limits>("xbjRange")) {}

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Proton inelastic form factors");
        desc.add<ParametersDescription>("structureFunctions", ParametersDescription().setName<int>(11));
        desc.add<ParametersDescription>("integrator", ParametersDescription().setName<std::string>("gsl"));
        desc.add<Limits>("xbjRange", Limits{1.e-9, 1.});
        return desc;
      }

    protected:
      void compute() override {
        setFEFM(integr_->integrate([this](double xbj) { return sf_->F2(xbj, q2_) * std::pow(xbj, -1); }, xbj_range_),
                integr_->integrate([this](double xbj) { return sf_->F2(xbj, q2_) * std::pow(xbj, -3); }, xbj_range_));
      }

    private:
      const std::unique_ptr<strfun::Parameterisation> sf_;
      const std::unique_ptr<AnalyticIntegrator> integr_;
      const Limits xbj_range_;
    };
  }  // namespace formfac
}  // namespace cepgen
using cepgen::formfac::InelasticFormFactors;
REGISTER_FORMFACTORS("InelasticNucleon", InelasticFormFactors);
