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

#include "CepGen/Event/Event.h"
#include "CepGen/Modules/SurvivalFactorFactory.h"
#include "CepGen/SurvivalFactors/SurvivalFactor.h"

using namespace std::string_literals;

namespace cepgen::survival {
  /// Constant survival factor, independent of kinematics
  class Constant final : public SurvivalFactor {
  public:
    explicit Constant(const ParametersList& params) : SurvivalFactor(params), value_(steer<double>("value")) {}

    static ParametersDescription description() {
      auto desc = SurvivalFactor::description();
      desc.setDescription("Constant survival factor");
      desc.add("value", 1.).setDescription("survival factor value");
      return desc;
    }

    double operator()(const Event&) const override { return value_; }

  private:
    const double value_;
  };
}  // namespace cepgen::survival
using cepgen::survival::Constant;
REGISTER_SURVIVAL_FACTOR("constant", Constant);
