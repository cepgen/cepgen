/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

#ifdef _WIN32
#include <libloaderapi.h>
#else
#include <dlfcn.h>
#endif

namespace cepgen::utils {
  std::atomic<int> gSignal;  ///< Abort signal handler
}

namespace cepgen {
  static const auto os_independent_path = [](const std::string& path, bool match) -> std::string {
#ifdef _WIN32
    return match ? path + ".dll" : path;
#elif defined(__APPLE__)
    return match ? "lib" + path + ".dylib" : path;
#else
    return match ? "lib" + path + ".so" : path;
#endif
  };

  bool loadLibrary(const std::string& path, bool match) {
    const auto full_path = os_independent_path(path, match);
    if (loaded_libraries.count(full_path) > 0)  // library is already loaded
      return true;
    if (callPath(full_path, [&full_path](const auto& file_path) -> bool {
#ifdef _WIN32
          if (auto* handle = LoadLibraryA(file_path.c_str()); handle != nullptr) {
            loaded_libraries.insert(std::make_pair(full_path, handle));
            return true;
          }
          CG_WARNING("loadLibrary") << "Failed to load library '" << file_path << "'.\n\t"
                                    << "Error code #" << GetLastError() << ".";
          return false;
#else
          if (auto* handle = ::dlopen(file_path.c_str(), RTLD_LAZY | RTLD_GLOBAL); handle != nullptr) {
            loaded_libraries.insert(std::make_pair(full_path, handle));
            return true;
          }
          const char* err = ::dlerror();
          CG_WARNING("loadLibrary") << "Failed to load library '" << file_path << "'."
                                    << (err != nullptr ? utils::format("\n\t%s", err) : "");
          return false;
#endif
        })) {
      CG_DEBUG("loadLibrary") << "Loaded library \"" << path << "\".";
      return true;
    }
    invalid_libraries.emplace_back(path);
    CG_DEBUG("loadLibrary") << "Library \"" << path << "\" (" << full_path << ") does not exist.";
    return false;
  }

  bool unloadLibrary(const std::string& path, bool match) {
    const auto full_path = os_independent_path(path, match);
    if (loaded_libraries.count(full_path) == 0) {  // library is not loaded
      CG_WARNING("unloadLibrary") << "Requested to remove library '" << full_path
                                  << "' from runtime environment, while it was not loaded. Libraries loaded:\n"
                                  << loaded_libraries << ".";
      return true;
    }
#ifdef _WIN32
    if (FreeLibrary(loaded_libraries.at(full_path))) {
      loaded_libraries.erase(full_path);
      CG_DEBUG("unloadLibrary") << "Successfully unloaded library '" << full_path
                                << "' from runtime environment. Remaining libraries loaded:\n"
                                << loaded_libraries << ".";
      return true;
    }
    CG_WARNING("unloadLibrary") << "Failed to unload library '" << full_path << "'.\n\t"
                                << "Error code #" << GetLastError() << ".";
#else
    if (const auto ret = ::dlclose(loaded_libraries.at(full_path)); ret == 0) {
      loaded_libraries.erase(full_path);
      CG_DEBUG("unloadLibrary") << "Successfully unloaded library '" << full_path
                                << "' from runtime environment. Remaining libraries loaded:\n"
                                << loaded_libraries << ".";
      return true;
    }
    const char* err = ::dlerror();
    CG_WARNING("unloadLibrary") << "Failed to unload library '" << full_path << "'."
                                << (err != nullptr ? utils::format("\n\t%s", err) : "");
#endif
    return false;
  }

  bool callPath(const std::string& local_path, const std::function<bool(const std::string&)>& callback) {
    if (search_paths.empty()) {
      CG_WARNING("callPath") << "List of search paths is empty.";
      return false;
    }
    for (const auto& search_path : search_paths)
      if (const auto the_path = fs::path{search_path} / local_path; utils::fileExists(the_path) && callback)
        return callback(the_path);
    return false;
  }

  void initialise(bool safe_mode) {
    // parse all particles properties
    static const std::string pdg_file = "";
    search_paths = utils::env::searchPaths();
    CG_DEBUG("initialise") << utils::s("Search path", search_paths.size(), false) << ": " << search_paths << ".";

    try {
      printHeader();  // display header message
    } catch (const Exception& e) {
      e.dump();
    }

    // particles table parsing
    if (!callPath("mass_width_2023.txt", [](const auto& path) {
          pdg::MCDFileParser::parse(path);
          return true;
        }))
      CG_WARNING("init") << "No particles definition file found.";
    if (PDG::get().size() < 10)
      CG_WARNING("init") << "Only " << utils::s("particle", PDG::get().size(), true)
                         << " are defined in the runtime environment.\n\t"
                         << "Make sure the path to the MCD file is correct.";

    std::string addons_file;
    for (const auto& search_path : search_paths) {
      if (const fs::path the_path{search_path}; addons_file.empty() && utils::fileExists(the_path / "CepGenAddOns.txt"))
        addons_file = the_path / "CepGenAddOns.txt";
      utils::env::append("LD_LIBRARY_PATH", search_path);
    }

    // load all necessary modules
    if (!safe_mode && !addons_file.empty())
      for (const auto& lib : utils::split(utils::readFile(addons_file), '\n'))
        loadLibrary(lib, true);
    loadLibrary("CepGenProcesses", true);
    if (!invalid_libraries.empty())
      CG_WARNING("init") << "Failed to load the following libraries:\n\t" << invalid_libraries << ".";

    // greeting message
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
    if (!callPath("README", [](const auto& path) {
          if (!utils::fileExists(path))
            return false;
          CG_LOG << utils::readFile(path);
          return true;
        }))
      CG_WARNING("printHeader") << "Failed to open README file.";
  }
}  // namespace cepgen
