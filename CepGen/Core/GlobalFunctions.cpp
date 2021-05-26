#include "CepGen/Generator.h"
#include "CepGen/Version.h"

#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Filesystem.h"

#include <fstream>
#include <atomic>

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
#else
    const auto fullpath = match ? "lib" + path + ".so" : path;
#endif
    if (!utils::fileExists(fullpath)) {
      CG_DEBUG("loadLibrary") << "Library \"" << path << "\" does not exist.";
      return false;
    }

#ifdef _WIN32
    if (LoadLibraryA(fullpath.c_str()) == nullptr) {
      CG_DEBUG("loadLibrary") << "Failed to load library \"" << path << "\".\n\t"
                              << "Error code #" << GetLastError() << ".";
      invalid_libraries.emplace_back(path);
      return false;
    }
#else
    if (dlopen(fullpath.c_str(), RTLD_LAZY | RTLD_GLOBAL) == nullptr) {
      const char* err = dlerror();
      CG_WARNING("loadLibrary") << "Failed to load library \"" << path << "\"."
                                << (err != nullptr ? utils::format("\n\t%s", err) : "");
      invalid_libraries.emplace_back(path);
      return false;
    }
#endif
    CG_DEBUG("loadLibrary") << "Loaded library \"" << path << "\".";
    loaded_libraries.emplace_back(path);
    return true;
  }

  void initialise(bool safe_mode) {
    //--- parse all particles properties
    static const std::string pdg_file = "";
    search_paths =
        std::vector<std::string>{utils::environ("CEPGEN_PATH", "."), fs::path() / "usr" / "share" / "CepGen"};

    //--- header message
    try {
      printHeader();
    } catch (const Exception& e) {
      e.dump();
    }

    //--- particles table parsing
    std::string mcd_file, addons_file;
    for (const auto& path : search_paths) {
      if (mcd_file.empty() && utils::fileExists(path + "/mass_width_2020.mcd"))
        mcd_file = path + "/mass_width_2020.mcd";
      if (addons_file.empty() && utils::fileExists(path + "/CepGenAddOns.txt"))
        addons_file = path + "/CepGenAddOns.txt";
    }
    if (mcd_file.empty())
      CG_WARNING("init") << "No particles definition file found.";
    else
      pdg::MCDFileParser::parse(mcd_file);
    if (PDG::get().size() < 10)
      CG_WARNING("init") << "Only " << utils::s("particle", PDG::get().size(), true)
                         << " are defined in the runtime environment.\n\t"
                         << "Make sure the path to the MCD file is correct.";

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
        CG_LOG("printHeader") << std::string(std::istreambuf_iterator<char>(hf), std::istreambuf_iterator<char>());
        return;
      }
    }
    throw CG_WARNING("printHeader") << "Failed to open README file.";
  }
}  // namespace cepgen
