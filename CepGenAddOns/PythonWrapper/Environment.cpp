/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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

// clang-format off
#include "CepGenAddOns/PythonWrapper/Environment.h"
#include "CepGenAddOns/PythonWrapper/Error.h"
// clang-format on

#include <algorithm>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace python {
    //------------------------------------------------------------------
    // Python API helpers
    //------------------------------------------------------------------

    Environment::Environment() {
      for (const auto& path : std::vector<std::string>{utils::env::get("CEPGEN_PATH", "."),
                                                       fs::current_path(),
                                                       fs::current_path() / "Cards",
                                                       fs::current_path().parent_path() / "Cards",
                                                       fs::current_path().parent_path().parent_path() / "Cards",
                                                       "/usr/share/CepGen/Cards"})
        utils::env::append("PYTHONPATH", path);
      CG_DEBUG("Python:Environment") << "PYTHONPATH set to " << utils::env::get("PYTHONPATH") << ".";

      Py_InitializeEx(1);
      if (!initialised())
        throw CG_FATAL("Python:Environment") << "Failed to initialise the Python environment!";
#if PY_VERSION_HEX >= 0x03080000
      PyConfig_InitPythonConfig(&config_);
#endif
      utils::env::set("PYTHONDONTWRITEBYTECODE", "1");
    }

    Environment::~Environment() {
      if (!initialised())
        CG_FATAL("Python:Environment")
            << "Python environment is set to be finalised while it was not initialised in the first place.";
      Py_Finalize();
    }

    bool Environment::initialised() { return Py_IsInitialized(); }

    void Environment::setProgramName(const std::string& filename) {
      const size_t fn_len = filename.length() + 1;
#ifdef PYTHON2
      char* sfilename = new char[fn_len];
      snprintf(sfilename, fn_len, "%s", filename.c_str());
      const std::string readable_s_filename(sfilename);
#else
      wchar_t* sfilename = new wchar_t[fn_len];
      swprintf(sfilename, fn_len, L"%s", filename.c_str());
      const std::wstring readable_s_filename(sfilename);
#endif
      if (!sfilename)
        throw CG_FATAL("PythonHandler") << "Invalid filename provided to the Python cards parser!";
#if PY_VERSION_HEX >= 0x03080000
      config_.program_name = sfilename;
#else
      Py_SetProgramName(sfilename);
#endif
      delete[] sfilename;
      CG_DEBUG("Python:setProgramName") << "Programme name set to \"" << readable_s_filename << "\".";
    }
  }  // namespace python
}  // namespace cepgen
