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
// clang-format on

#include <algorithm>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Filesystem.h"

namespace cepgen {
  namespace python {
    //------------------------------------------------------------------
    // Python API helpers
    //------------------------------------------------------------------

    Environment::Environment(const ParametersList& params) : SteeredObject(params) {
      const auto cepgen_path = fs::path(utils::env::get("CEPGEN_PATH", "."));
      for (const auto& path : std::vector<std::string>{cepgen_path,
                                                       cepgen_path / "python",
                                                       cepgen_path / "python_modules",
                                                       cepgen_path / "build" / "python",
                                                       cepgen_path / "build" / "python_modules",
                                                       fs::current_path(),
                                                       fs::current_path() / "python",
                                                       fs::current_path() / "python_modules",
                                                       "/usr/share/CepGen/python",
                                                       "/usr/share/CepGen/python_modules"})
        utils::env::append("PYTHONPATH", path);
      CG_DEBUG("Python:Environment") << "PYTHONPATH set to " << utils::env::get("PYTHONPATH") << ".";

#if PY_VERSION_HEX >= 0x03080000
      PyConfig_InitPythonConfig(&config_);
      config_.parser_debug = steer<int>("debug");
      config_.verbose = steer<int>("verbosity");
      Py_InitializeFromConfig(&config_);
#else
      Py_DebugFlag = steer<int>("debug");
      Py_VerboseFlag = steer<int>("verbosity");
      Py_InitializeEx(1);
#endif
      if (!initialised())
        throw CG_FATAL("Python:Environment") << "Failed to initialise the Python environment!";
      utils::env::set("PYTHONDONTWRITEBYTECODE", "1");
      if (const auto& name = steer<std::string>("name"); !name.empty())
        setProgramName(name);
    }

    Environment::~Environment() {
      if (!initialised())
        CG_WARNING("Python:Environment")
            << "Python environment is set to be finalised while it was not initialised in the first place.";
      else
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

    ParametersDescription Environment::description() {
      auto desc = ParametersDescription();
      desc.add<int>("verbosity", 0).setDescription("overall Python verbosity");
      desc.add<int>("debug", 0).setDescription("debugging level");
      return desc;
    }
  }  // namespace python
}  // namespace cepgen
