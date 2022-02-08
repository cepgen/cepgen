/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include <matheval.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    class FunctionalLibmatheval final : public Functional {
    public:
      explicit FunctionalLibmatheval(const ParametersList&);

      static ParametersDescription description();

      double eval(const std::vector<double>&) const;

    private:
      struct eval_deleter {
        void operator()(void* eval) { evaluator_destroy(eval); }
      };
      std::unique_ptr<void, eval_deleter> eval_;
      std::vector<std::string> parsed_vars_;
      char** c_parsed_vars_;
    };

    FunctionalLibmatheval::FunctionalLibmatheval(const ParametersList& params) : Functional(params) {
      eval_.reset(evaluator_create(const_cast<char*>(expression_.c_str())));
      if (!eval_)
        throw CG_ERROR("FunctionalLibmatheval") << "Evaluator was not properly initialised!";
      int num_vars;
      evaluator_get_variables(eval_.get(), &c_parsed_vars_, &num_vars);
      for (int i = 0; i < num_vars; ++i)
        parsed_vars_.emplace_back(c_parsed_vars_[i]);
      if (parsed_vars_.size() != vars_.size())
        throw CG_FATAL("FunctionalLibmatheval")
            << "Parsed " << utils::s("variable", num_vars, true) << ": " << parsed_vars_ << " where "
            << utils::s("variable", vars_.size(), true) << " is/are expected: " << vars_ << "!";
    }

    double FunctionalLibmatheval::eval(const std::vector<double>& x) const {
      if (parsed_vars_.size() != x.size())
        throw CG_FATAL("FunctionalLibmatheval") << "Invalid number of variables fed to the evaluator!";
      return evaluator_evaluate(eval_.get(), parsed_vars_.size(), c_parsed_vars_, const_cast<double*>(x.data()));
    }

    ParametersDescription FunctionalLibmatheval::description() {
      auto desc = Functional::description();
      desc.setDescription("libmatheval evaluator");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_FUNCTIONAL("libmatheval", FunctionalLibmatheval)
