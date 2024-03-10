/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Caller.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    Caller::Caller() : oldcout_(std::cout.rdbuf(os_cout_.rdbuf())), oldcerr_(std::cerr.rdbuf(os_cerr_.rdbuf())) {}

    Caller::~Caller() {
      std::cout.rdbuf(oldcout_);
      std::cerr.rdbuf(oldcerr_);
      if (const auto& str = os_cout_.str(); !str.empty())
        CG_DEBUG("Caller") << "At end of caller call, the following output was generated:\n" << str;
      if (const auto& str = os_cerr_.str(); !str.empty())
        CG_WARNING("Caller") << "At end of caller call, the following error stream was generated:\n" << str;
    }

    std::string Caller::output() const { return os_cout_.str(); }

    std::string Caller::error() const { return os_cerr_.str(); }

    std::string Caller::call(const std::vector<std::string>& commands) { return call(utils::merge(commands, " ")); }

    std::string Caller::call(const std::string& command) {
      auto pipe = popen(command.c_str(), "r");
      if (!pipe)
        throw CG_FATAL("Caller") << "Failed to call the command '" << command << "'.";
      std::array<char, 128> buffer;
      std::string out;
      while (!feof(pipe))
        if (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
          std::cout << buffer.data();
          out += buffer.data();
        }
      auto rc = pclose(pipe);
      if (rc != EXIT_SUCCESS)
        throw CG_FATAL("Caller") << "Command '" << command << "' failed with return code '" << rc << "'.";
      return out;
    }
  }  // namespace utils
}  // namespace cepgen
