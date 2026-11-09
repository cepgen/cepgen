/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2026  Laurent Forthomme
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

#include <CepGen/Core/Exception.h>
#include <CepGen/Core/ParametersList.h>
#include <pybind11/stl.h>

#include "CepGenPython/PythonUtils.h"

namespace py = pybind11;

namespace {
  std::string type_name(const py::handle& object) {
    return object.attr("__class__").attr("__name__").cast<std::string>();
  }
}  // namespace

namespace cepgen::python {
  ParametersList py_dict_to_plist(const py::dict& dict) {
    ParametersList plist;
    for (const auto& item : dict) {
      const std::string key = py::str(item.first);
      const auto value = item.second;
      const auto value_type = type_name(value);
      if (value_type == "int")
        plist.set(key, value.cast<int>());
      else if (value_type == "str")
        plist.set(key, value.cast<std::string>());
      else if (value_type == "float")
        plist.set(key, value.cast<double>());
      else if (value_type == "dict")
        plist.set(key, py_dict_to_plist(value.cast<py::dict>()));
      else if (value_type == "tuple") {
        const auto tuple = value.cast<py::tuple>();
        const auto element_type = type_name(tuple[0]);
        if (element_type == "float") {
          if (const auto items = value.cast<std::vector<double> >(); items.size() == 2)
            plist.set(key, Limits{items.at(0), items.at(1)});
        } else
          throw CG_FATAL("py_dict_to_plist") << "Tuple unpacking is not (yet) handling the Python '" << element_type
                                             << "' type for key='" << key << "'.";
      } else if (value_type == "list") {
        // find the type of list from its first element
        const auto list = value.cast<py::list>();
        const auto element_type = type_name(list[0]);
        if (element_type == "int")
          plist.set(key, value.cast<std::vector<int> >());
        else if (element_type == "str")
          plist.set(key, value.cast<std::vector<std::string> >());
        else if (element_type == "float")
          plist.set(key, value.cast<std::vector<double> >());
        else if (element_type == "dict") {
          std::vector<ParametersList> plists;
          for (const auto& element : list)
            plists.emplace_back(py_dict_to_plist(element.cast<py::dict>()));
          plist.set(key, plists);
        } else
          throw CG_FATAL("py_dict_to_plist")
              << "Failed to unpack a Python list for elements of '" << element_type << "' type.";
      } else
        throw CG_FATAL("py_dict_to_plist")
            << "Failed to unpack a Python '" << value_type << "' type for key='" << key << "'.";
    }
    return plist;
  }

  py::dict plist_to_py_dict(const ParametersList& plist) {
    py::dict dict;
    // a key may be stored under several types: the first type handled below takes precedence
    const auto set_default = dict.attr("setdefault");
    for (const auto& key : plist.keysOf<int>())
      set_default(key, plist.get<int>(key));
    for (const auto& key : plist.keysOf<std::string>())
      set_default(key, plist.get<std::string>(key));
    for (const auto& key : plist.keysOf<double>())
      set_default(key, plist.get<double>(key));
    for (const auto& key : plist.keysOf<Limits>()) {
      const auto& limits = plist.get<Limits>(key);
      set_default(key, py::make_tuple(limits.min(), limits.max()));
    }
    for (const auto& key : plist.keysOf<ParametersList>())
      set_default(key, plist_to_py_dict(plist.get<ParametersList>(key)));
    for (const auto& key : plist.keysOf<std::vector<int> >())
      set_default(key, plist.get<std::vector<int> >(key));
    for (const auto& key : plist.keysOf<std::vector<std::string> >())
      set_default(key, plist.get<std::vector<std::string> >(key));
    for (const auto& key : plist.keysOf<std::vector<double> >())
      set_default(key, plist.get<std::vector<double> >(key));
    return dict;
  }
}  // namespace cepgen::python
