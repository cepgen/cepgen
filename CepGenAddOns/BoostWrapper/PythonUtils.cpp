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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGenAddOns/BoostWrapper/PythonUtils.h"

cepgen::ParametersList py_dict_to_plist(const py::dict& dict) {
  cepgen::ParametersList plist;
  for (auto it = py::stl_input_iterator<py::tuple>(dict.attr("items")()); it != py::stl_input_iterator<py::tuple>();
       ++it) {
    py::tuple kv = *it;
    const std::string key{py::extract<const char*>(py::str(kv[0]))};
    py::extract<py::object> val_ext(kv[1]);
    const std::string val_type{py::extract<const char*>(val_ext().attr("__class__").attr("__name__"))};
    if (val_type == "int")
      plist.set<int>(key, py::extract<int>(val_ext()));
    else if (val_type == "str")
      plist.set<std::string>(key, std::string{py::extract<const char*>(val_ext())});
    else if (val_type == "float")
      plist.set<double>(key, py::extract<double>(val_ext()));
    else if (val_type == "dict")
      plist.set<cepgen::ParametersList>(key, py_dict_to_plist(py::extract<py::dict>(val_ext())));
    else if (val_type == "list") {
      // find the type of list from its first element
      const py::list& list = py::extract<py::list>(val_ext());
      const std::string el_type{py::extract<const char*>(list[0].attr("__class__").attr("__name__"))};
      if (el_type == "int")
        plist.set(key, py_list_to_std_vector<int>(py::extract<py::list>(val_ext())));
      else if (el_type == "str")
        plist.set(key, py_list_to_std_vector<std::string>(py::extract<py::list>(val_ext())));
      else if (el_type == "float")
        plist.set(key, py_list_to_std_vector<double>(py::extract<py::list>(val_ext())));
      else
        throw CG_FATAL("py_dict_to_plist")
            << "Failed to unpack a Python list for elements of '" << val_type << "' type.";
    } else
      throw CG_FATAL("py_dict_to_plist") << "Failed to unpack a Python '" << val_type << "' type for key='" << key
                                         << "'.";
  }
  return plist;
}
