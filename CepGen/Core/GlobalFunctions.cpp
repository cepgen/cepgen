/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#include <atomic>
#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

#ifdef _WIN32
#include <libloaderapi.h>
#else
#include <dlfcn.h>
#endif

namespace cepgen {
  namespace utils {
    std::atomic<int> gSignal;  ///< Abort signal handler
  }                            // namespace utils

  bool loadLibrary(const std::string& path, bool match) {
    if (utils::contains(loaded_libraries, path))
      return true;
#ifdef _WIN32
    const auto fullpath = match ? path + ".dll" : path;
#elif defined(__APPLE__)
    const auto fullpath = match ? "lib" + path + ".dylib" : path;
#else
    const auto fullpath = match ? "lib" + path + ".so" : path;
#endif
    for (const auto& search_path : search_paths) {
      fs::path the_path{search_path};
      the_path /= fullpath;
      if (!utils::fileExists(the_path))
        continue;

#ifdef _WIN32
      if (LoadLibraryA(the_path.c_str()) == nullptr) {
        CG_DEBUG("loadLibrary") << "Failed to load library \"" << the_path << "\".\n\t"
                                << "Error code #" << GetLastError() << ".";
        invalid_libraries.emplace_back(path);
        return false;
      }
#else
      if (dlopen(the_path.c_str(), RTLD_LAZY | RTLD_GLOBAL) == nullptr) {
        const char* err = dlerror();
        CG_WARNING("loadLibrary") << "Failed to load library " << the_path << "."
                                  << (err != nullptr ? utils::format("\n\t%s", err) : "");
        invalid_libraries.emplace_back(path);
        return false;
      }
#endif
      CG_DEBUG("loadLibrary") << "Loaded library \"" << path << "\".";
      loaded_libraries.emplace_back(path);
      return true;
    }
    CG_DEBUG("loadLibrary") << "Library \"" << path << "\" does not exist.";
    return false;
  }

  void initialise(bool safe_mode) {
    //--- parse all particles properties
    static const std::string pdg_file = "";
    search_paths = utils::env::searchPaths();

    //--- particles table parsing
    std::string mcd_file, addons_file;
    for (const auto& path : search_paths) {
      if (mcd_file.empty() && utils::fileExists(path + "/mass_width_2021.mcd"))
        mcd_file = path + "/mass_width_2021.mcd";
      if (addons_file.empty() && utils::fileExists(path + "/CepGenAddOns.txt"))
        addons_file = path + "/CepGenAddOns.txt";
      utils::env::append("LD_LIBRARY_PATH", path);
    }
    if (mcd_file.empty())
      CG_WARNING("init") << "No particles definition file found.";
    else
      pdg::MCDFileParser::parse(mcd_file);
    if (PDG::get().size() < 10)
      CG_WARNING("init") << "Only " << utils::s("particle", PDG::get().size(), true)
                         << " are defined in the runtime environment.\n\t"
                         << "Make sure the path to the MCD file is correct.";

    //--- header message
    try {
      printHeader();
    } catch (const Exception& e) {
      e.dump();
    }

    //--- load all necessary modules
    if (!safe_mode && !addons_file.empty()) {
      std::ifstream addons(addons_file);
      std::string lib;
      while (std::getline(addons, lib))
        loadLibrary(lib, true);
    }
    loadLibrary("CepGenProcesses", true);
    if (!invalid_libraries.empty())
      CG_WARNING("init") << "Failed to load the following libraries:\n\t" << invalid_libraries << ".";

    //--- greeting message
    CG_INFO("init").log([&](auto& log) {
      log << "CepGen " << version::tag << " (" << version::extended << ") "
          << "initialised";
      if (!loaded_libraries.empty())
        log << " with " << utils::s("add-on", loaded_libraries.size(), true) << ":\n\t" << loaded_libraries << ".\n\t";
      else
        log << ". ";
      log << "Greetings!";
    });
  }

  void printHeader() {
    for (const auto& path : search_paths) {
      std::ifstream hf(path + "/README");
      if (hf.good()) {
        CG_LOG << std::string(std::istreambuf_iterator<char>(hf), std::istreambuf_iterator<char>());
        return;
      }
    }
    CG_WARNING("printHeader") << "Failed to open README file.";
  }
}  // namespace cepgen
