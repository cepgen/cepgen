/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/SurvivalFactorFactory.h"
#include "CepGen/SurvivalFactors/SurvivalFactor.h"
#include "CepGen/Utils/Functional.h"

using namespace std::string_literals;

namespace cepgen::survival {
  /// Survival factor parameterised as a mathematical expression of event-level quantities
  /// \note Variables follow the event browser syntax, e.g. m(4) for the two-parton system mass, pt(41) for the
  ///  positive-z parton transverse momentum, or m(3) for the positive-z outgoing beam system mass.
  class FunctionalSurvivalFactor final : public SurvivalFactor {
  public:
    explicit FunctionalSurvivalFactor(const ParametersList& params)
        : SurvivalFactor(params), variables_(steer<std::vector<std::string> >("variables")) {
      if (variables_.empty())
        throw CG_FATAL("SurvivalFactor:functional") << "At least one event variable must be specified.";
      functional_ = FunctionalFactory::get().build(
          steer<std::string>("functionalEngine"),
          utils::Functional::fromExpression(steer<std::string>("expression"), variables_));
      values_.resize(variables_.size());
    }

    static ParametersDescription description() {
      auto desc = SurvivalFactor::description();
      desc.setDescription("Functional survival factor");
      desc.add("expression", "1."s).setDescription("mathematical expression of the survival factor");
      desc.add("variables", std::vector<std::string>{"m(4)"})
          .setDescription("event variables (event browser syntax) used in the expression");
      desc.add("functionalEngine", "root"s).setDescription("functional parser to use");
      return desc;
    }

    double operator()(const Event& event) const override {
      for (size_t i = 0; i < variables_.size(); ++i)
        values_[i] = browser_.get(event, variables_[i]);
      return (*functional_)(values_);
    }

  private:
    const std::vector<std::string> variables_;
    const utils::EventBrowser browser_;
    std::unique_ptr<utils::Functional> functional_;
    mutable std::vector<double> values_;
  };
}  // namespace cepgen::survival
using cepgen::survival::FunctionalSurvivalFactor;
REGISTER_SURVIVAL_FACTOR("functional", FunctionalSurvivalFactor);
