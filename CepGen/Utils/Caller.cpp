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

#include <memory>
#include <string>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Caller.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    std::string Caller::call(const std::vector<std::string>& commands) { return call(utils::merge(commands, " ")); }

    std::string Caller::call(const std::string& command) {
      auto pipe = popen(command.c_str(), "r");
      if (!pipe)
        throw CG_FATAL("Caller") << "Failed to call the command '" << command << "'.";
      std::array<char, 128> buffer;
      std::string result;
      while (!feof(pipe))
        if (fgets(buffer.data(), buffer.size(), pipe) != nullptr)
          result += buffer.data();
      auto rc = pclose(pipe);
      CG_DEBUG("Caller").log([&](auto& log) {
        log << "Command '" << command << "' returned code '" << rc << "'.";
        if (!result.empty())
          log << " Message: '" << result << "'.";
      });
      if (rc != EXIT_SUCCESS)
        throw CG_FATAL("Caller") << "Command '" << command << "' failed with return code '" << rc << "'.";
      return result;
    }
  }  // namespace utils
}  // namespace cepgen
