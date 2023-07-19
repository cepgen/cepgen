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

#include "CepGen/Generator.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGenAddOns/BoostWrapper/PythonUtils.h"

namespace {
  BOOST_PYTHON_MODULE(pycepgen) {
    cepgen::initialise();

    py::class_<cepgen::strfun::Parameterisation, boost::noncopyable>(
        "_StructureFunctions", "nucleon structure functions modelling", py::no_init)
        .add_static_property("name",
                             py::make_function(&cepgen::strfun::Parameterisation::name,
                                               py::return_value_policy<py::copy_const_reference>()))
        .def("F2", &cepgen::strfun::Parameterisation::F2)
        .def("FL", &cepgen::strfun::Parameterisation::FL)
        .def("F1", &cepgen::strfun::Parameterisation::F1);

    EXPOSE_FACTORY(cepgen::StructureFunctionsFactory,
                   int,
                   "StructureFunctionsFactory",
                   "a structure functions evaluator objects factory");
  }
}  // namespace
