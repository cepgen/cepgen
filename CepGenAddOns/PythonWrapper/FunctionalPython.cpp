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

#include <Python.h>

#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"

namespace cepgen {
  namespace utils {
    class FunctionalPython final : public Functional {
    public:
      explicit FunctionalPython(const ParametersList&);

      static ParametersDescription description();

      double eval(const std::vector<double>&) const;

    private:
    };

    FunctionalPython::FunctionalPython(const ParametersList& params) : Functional(params) {}

    double FunctionalPython::eval(const std::vector<double>& x) const { return 0.; }

    ParametersDescription FunctionalPython::description() {
      auto desc = Functional::description();
      desc.setDescription("Python mathematical expression evaluator");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_FUNCTIONAL("python", FunctionalPython)
