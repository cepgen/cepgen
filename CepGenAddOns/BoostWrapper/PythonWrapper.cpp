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

#include <boost/mpl/vector.hpp>
#include <boost/python.hpp>

#include "CepGen/Generator.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

#define EXPOSE_FACTORY(obj, key, name, description)                                                             \
  py::class_<obj, boost::noncopyable>(name, description, py::no_init)                                           \
      .def("build", adapt_unique(+[](const key& mod) { return obj::get().build(mod); }))                        \
      .def("build", adapt_unique(+[](const cepgen::ParametersList& plist) { return obj::get().build(plist); })) \
      .add_static_property(                                                                                     \
          "modules", +[]() { return std_vector_to_py_list(obj::get().modules()); })

namespace {
  namespace py = boost::python;

  template <class T>
  py::list std_vector_to_py_list(const std::vector<T>& vec) {
    py::list lst;
    std::for_each(vec.begin(), vec.end(), [&](const auto& t) { lst.append(t); });
    return lst;
  }

  template <typename T, typename... Args>
  py::object adapt_unique(std::unique_ptr<T> (*fn)(Args...)) {
    return py::make_function([fn](Args... args) { return fn(args...).release(); },
                             py::return_value_policy<py::manage_new_object>(),
                             boost::mpl::vector<T*, Args...>());
  }

  template <typename T, typename C, typename... Args>
  py::object adapt_unique(std::unique_ptr<T> (C::*fn)(Args...)) {
    return py::make_function([fn](C& self, Args... args) { return (self.*fn)(args...).release(); },
                             py::return_value_policy<py::manage_new_object>(),
                             boost::mpl::vector<T*, C&, Args...>());
  }

  BOOST_PYTHON_MODULE(pycepgen) {
    cepgen::initialise();

    py::class_<cepgen::strfun::Parameterisation, boost::noncopyable>(
        "_StructureFunctions", "nucleon structure functions modelling", py::no_init)
        .def("F2", &cepgen::strfun::Parameterisation::F2)
        .def("FL", &cepgen::strfun::Parameterisation::FL)
        .def("F1", &cepgen::strfun::Parameterisation::F1);

    EXPOSE_FACTORY(cepgen::StructureFunctionsFactory,
                   int,
                   "StructureFunctionsFactory",
                   "a structure functions evaluator objects factory");
  }
}  // namespace
