/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include <apfel/apfelxx.h>

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"

namespace cepgen::apfelpp {
  class AlphaS final : public cepgen::Coupling {
  public:
    explicit AlphaS(const ParametersList& params)
        : cepgen::Coupling(params),
          use_tabulated_(steer<bool>("useTabulated")),
          alpha_s_(new apfel::AlphaQCD(steer<double>("alphaSref"),
                                       steer<double>("muQCDref"),
                                       steer<std::vector<double> >("quarkThresholds"),
                                       steer<int>("order"))),
          tab_params_(steer<ParametersList>("tabulatedParameters")),
          alpha_s_tab_(use_tabulated_ ? new apfel::TabulateObject<double>(*alpha_s_,
                                                                          tab_params_.get<int>("numValues"),
                                                                          tab_params_.get<Limits>("Qrange").min(),
                                                                          tab_params_.get<Limits>("Qrange").max(),
                                                                          tab_params_.get<int>("order"),
                                                                          tab_params_.get<double>("Lambda"))
                                      : nullptr) {
      apfel::Banner();
    }

    static ParametersDescription description() {
      auto desc = cepgen::Coupling::description();
      desc.setDescription("APFEL++ alpha(S) evolution algorithm");
      desc.add<bool>("useTabulated", true).setDescription("use the tabulated, fast values interpolator?");
      desc.add<double>("alphaSref", 0.118);
      desc.add<double>("muQCDref", 91.1876);
      desc.add<std::vector<double> >("quarkThresholds", {0., 0., 0., M_SQRT2, 4.5, 175.});
      desc.add<int>("order", 2)
          .setDescription("QCD perturbative evolution order")
          .allow(0, "LO")
          .allow(1, "NLO")
          .allow(2, "NNLO")
          .allow(3, "NNNLO");

      auto tab_desc = ParametersDescription();
      tab_desc.add<int>("numValues", 100).setDescription("number of values evaluated to build the interpolation");
      tab_desc.add<Limits>("Qrange", {0.9, 1001.}).setDescription("Q range for the interpolation");
      tab_desc.add<int>("order", 3).setDescription("interpolation order");
      tab_desc.add<double>("Lambda", 0.25)
          .setDescription("Lambda parameter in the tabulation function (ln(ln(Q^2/Lambda^2))");
      desc.add("tabulatedParameters", tab_desc);
      return desc;
    }

    inline double operator()(double q) const override {
      return use_tabulated_ ? alpha_s_tab_->Evaluate(q) : alpha_s_->Evaluate(q);
    }

  private:
    const bool use_tabulated_;
    const std::unique_ptr<apfel::AlphaQCD> alpha_s_;
    const ParametersList tab_params_;
    const std::unique_ptr<apfel::TabulateObject<double> > alpha_s_tab_;
  };
}  // namespace cepgen::apfelpp
using AlphaSAPFELpp = cepgen::apfelpp::AlphaS;
REGISTER_ALPHAS_MODULE("apfelpp", AlphaSAPFELpp);
