/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025  Laurent Forthomme
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

#include <atmsp.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

namespace cepgen::atmsp {
  class Functional final : public cepgen::utils::Functional {
  public:
    explicit Functional(const ParametersList& params)
        : cepgen::utils::Functional(params), byte_code_(new ATMSB<double>) {
      ATMSP<double> parser;  // parsing/bytecode generation with error check
      if (const auto err = parser.parse(*byte_code_, expression_, utils::merge(vars_, ", ")); err)
        throw CG_ERROR("atmsp:Functional") << "Evaluator was not properly initialised. ATMSP error:\n"
                                           << parser.errMessage(err);
    }
    double eval() const override {
      for (size_t i = 0; i < values_.size(); ++i)  // set variable values
        byte_code_->var[i] = values_.at(i);
      return byte_code_->run();
    }

    static ParametersDescription description() {
      auto desc = cepgen::utils::Functional::description();
      desc.setDescription("ATMSP evaluator");
      return desc;
    }

  private:
    const std::unique_ptr<ATMSB<double> > byte_code_;  ///< Bytecode instance with SAME basic type as the parser
  };
}  // namespace cepgen::atmsp
using ATMSPFunctional = cepgen::atmsp::Functional;
REGISTER_FUNCTIONAL("atmsp", ATMSPFunctional);
