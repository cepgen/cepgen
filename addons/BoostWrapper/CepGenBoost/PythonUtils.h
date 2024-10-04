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

#ifndef CepGenBoost_PythonUtils_h
#define CepGenBoost_PythonUtils_h

#include <boost/mpl/vector.hpp>
#include <boost/python.hpp>

#define EXPOSE_FACTORY(obj, key, name, description)                                                               \
  py::class_<obj, boost::noncopyable>(name, description, py::no_init)                                             \
      .def("build", adapt_unique(+[](const key& mod) { return obj::get().build(mod); }))                          \
      .def("build", adapt_unique(+[](const key& mod, const py::dict& dict) {                                      \
             return obj::get().build(mod, py_dict_to_plist(dict));                                                \
           }))                                                                                                    \
      .def("build", adapt_unique(+[](const py::dict& dict) { return obj::get().build(py_dict_to_plist(dict)); })) \
      .def(                                                                                                       \
          "describe",                                                                                             \
          +[](const key& mod) {                                                                                   \
            std::ostringstream os;                                                                                \
            os << obj::get().describeParameters(mod);                                                             \
            return os.str();                                                                                      \
          })                                                                                                      \
      .add_static_property("modules", +[]() { return std_vector_to_py_list(obj::get().modules()); })

namespace cepgen {
  class ParametersList;
}

namespace py = boost::python;

template <class T>
py::list std_vector_to_py_list(const std::vector<T>& vec) {
  py::list list;
  std::for_each(vec.begin(), vec.end(), [&list](const auto& t) { list.append(t); });
  return list;
}

template <class T>
std::vector<T> py_tuple_to_std_vector(const py::tuple& tuple) {
  std::vector<T> vec;
  for (ssize_t i = 0; i < py::len(tuple); ++i)
    vec.emplace_back(py::extract<T>(tuple[i]));
  return vec;
}

template <class T>
std::vector<T> py_list_to_std_vector(const py::list& list) {
  std::vector<T> vec;
  for (ssize_t i = 0; i < py::len(list); ++i)
    vec.emplace_back(py::extract<T>(list[i]));
  return vec;
}

template <class T>
py::tuple std_vector_to_py_tuple(const std::vector<T>& vec) {
  return py::tuple(std_vector_to_py_list(vec));
}

cepgen::ParametersList py_dict_to_plist(const py::dict&);
py::dict plist_to_py_dict(const cepgen::ParametersList&);

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

template <typename T>
inline py::object adapt_reference(T* ptr) {
  typename py::reference_existing_object::apply<T*>::type converter;
  return py::object(py::handle(converter(ptr)));
}

#endif
