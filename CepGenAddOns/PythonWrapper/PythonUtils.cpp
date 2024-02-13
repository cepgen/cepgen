/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

namespace cepgen {
  namespace python {
    std::string pythonPath(const std::string& file) {
      const auto dir = fs::path{file}.remove_filename();
      if (!dir.empty()) {
        CG_DEBUG("Python") << "Adding {" << dir << "} to the default search paths.";
        utils::env::append("PYTHONPATH", dir);
      }

      const auto filename = utils::replace_all(fs::path{file}.replace_extension("").string() /* remove the extension */,
                                               {{"../", ".."}, {"/", "."}});
      CG_DEBUG("Python") << "Python path: " << filename;
      return filename;
    }

    std::vector<std::wstring> info() {
      auto* py_home = Py_GetPythonHome();
#ifdef PYTHON2
      std::wstring path{utils::towstring(std::string(Py_GetPath()))},
          home{utils::towstring(std::string(py_home ? py_home : "(not set)"))};
#else
      std::wstring path{Py_GetPath()}, home{py_home ? py_home : L"(not set)"};
#endif
      return std::vector<std::wstring>{
          utils::towstring("Python version: " + utils::replace_all(std::string{Py_GetVersion()}, "\n", " ")),
          utils::towstring("Platform: " + std::string(Py_GetPlatform())),
          utils::towstring("Home directory: ") + home,
          utils::towstring("Parsed path: ") + path};
    }
  }  // namespace python
}  // namespace cepgen
