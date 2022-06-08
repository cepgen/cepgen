/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    namespace env {
      std::string get(const std::string& var, const std::string& def) {
        const auto out = std::getenv(var.c_str());
        if (!out)
          return def;
        return std::string(out);
      }

      std::vector<std::string> searchPaths() {
        return std::vector<std::string>{env::get("CEPGEN_PATH", "."),
                                        fs::path() / "/usr" / "share" / "CepGen",
                                        fs::current_path().parent_path(),
                                        fs::current_path().parent_path().parent_path()};
      }

      void set(const std::string& var, const std::string& value) { setenv(var.c_str(), value.c_str(), 1); }

#ifdef _WIN32
      static constexpr char PATH_DELIM = ';';
#else
      static constexpr char PATH_DELIM = ':';
#endif

      void append(const std::string& var, const std::string& value) {
        auto env = split(env::get(var, ""), PATH_DELIM);
        env.emplace_back(value);
        normalise(env);
        setenv(var.c_str(), merge(env, std::string(1, PATH_DELIM)).c_str(), 1);
      }

      void unset(const std::string& var) { unsetenv(var.c_str()); }
    }  // namespace env
  }    // namespace utils
}  // namespace cepgen
