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
#include "CepGenPython/Environment.h"
#include "CepGenPython/Utils.h"

namespace cepgen::python {
  std::string pythonPath(const std::string& file) {
    const auto dir = fs::path{file}.remove_filename();
    if (!dir.empty()) {
      CG_DEBUG("Python") << "Adding {" << dir << "} to the default search paths.";
      utils::env::append("PYTHONPATH", dir);
    }

    const auto filename = utils::replaceAll(fs::path{file}.replace_extension("").string() /* remove the extension */,
                                            {{"../", ".."}, {"/", "."}});
    CG_DEBUG("Python") << "Python path: " << filename;
    return filename;
  }

  std::vector<std::wstring> info() {
#ifdef PYTHON2
    auto* py_home = Py_GetPythonHome();
    std::wstring path{utils::toWstring(std::string(Py_GetPath()))},
        home{utils::toWstring(std::string(py_home ? py_home : "(not set)"))};
#elif PY_VERSION_HEX < 0x030d0000  // python < 3.13
    auto* py_home = Py_GetPythonHome();
    std::wstring path{Py_GetPath()}, home{py_home ? py_home : L"(not set)"};
#else
    Environment env{ParametersList()};
    std::wstring path;
    for (Py_ssize_t i = 0; i < env.configuration().module_search_paths.length; ++i)
      path += std::wstring{env.configuration().module_search_paths.items[i]};
    std::wstring home{env.configuration().home};
    throw CG_ERROR("");
#endif
    return std::vector<std::wstring>{
        utils::toWstring("Python version: " + utils::replaceAll(std::string{Py_GetVersion()}, "\n", " ")),
        utils::toWstring("Platform: " + std::string(Py_GetPlatform())),
        utils::toWstring("Home directory: ") + home,
        utils::toWstring("Parsed path: ") + path};
  }
}  // namespace cepgen::python
