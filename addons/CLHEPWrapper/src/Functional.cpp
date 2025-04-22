/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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

#include <CLHEP/Evaluator/Evaluator.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"

using namespace std::string_literals;

namespace cepgen::clhep {
  class Functional final : public utils::Functional {
  public:
    explicit Functional(const ParametersList& params) : utils::Functional(params), evaluator_(new HepTool::Evaluator) {
      if (steer<bool>("useStdMath"))
        evaluator_->setStdMath();
    }

    static ParametersDescription description() {
      auto desc = utils::Functional::description();
      desc.setDescription("CLHEP functional evaluator");
      desc.add("useStdMath", true).setDescription("use the STL math environment?");
      return desc;
    }

    inline double eval() const override {
      for (size_t i = 0; i < vars_.size(); ++i)
        evaluator_->setVariable(vars_.at(i).data(), values_.at(i));
      const auto res = evaluator_->evaluate(expression_.data());
      if (evaluator_->status() != HepTool::Evaluator::OK)
        throw CG_ERROR("clhep:Functional")
            << "Error encountered while evaluating the expression:\n"
            << "  " << expression_ << "\n  "
            << (evaluator_->error_position() > 0 ? std::string(evaluator_->error_position(), '-') : ""s) << "^\n"
            << evaluator_->error_name();
      return res;
    }

  private:
    const std::unique_ptr<HepTool::Evaluator> evaluator_;
  };
}  // namespace cepgen::clhep
using CLHEPFunctional = cepgen::clhep::Functional;
REGISTER_FUNCTIONAL("clhep", CLHEPFunctional);
