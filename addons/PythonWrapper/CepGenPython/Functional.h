/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

#ifndef CepGenPython_Functional_h
#define CepGenPython_Functional_h

#include "CepGen/Utils/Functional.h"
#include "CepGenPython/Environment.h"
#include "CepGenPython/ObjectPtr.h"

namespace cepgen::python {
  class Functional final : public utils::Functional {
  public:
    explicit Functional(const ParametersList&);
    explicit Functional(const ObjectPtr&);

    double eval() const override;

    static ParametersDescription description();

  private:
    const std::unique_ptr<Environment> environment_;
    ObjectPtr mod_{nullptr};
    ObjectPtr func_{nullptr};
  };
}  // namespace cepgen::python

#endif
