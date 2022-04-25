/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#ifndef CepGen_Utils_Piper_h
#define CepGen_Utils_Piper_h

#include <memory>
#include <string>
#include <vector>

namespace cepgen {
  namespace utils {
    class Piper {
    public:
      explicit Piper(const std::string& command) : file_(popen(command.c_str(), "w"), pclose) {}

      struct Commands : std::vector<std::string> {
        using std::vector<std::string>::vector;
        Commands& operator+=(const std::vector<std::string>& oth) {
          std::copy(oth.begin(), oth.end(), std::back_inserter(*this));
          return *this;
        }
        Commands& operator+=(const std::string& str) {
          emplace_back(str);
          return *this;
        }
        friend std::ostream& operator<<(std::ostream& os, const Commands& cmds) {
          std::string sep;
          os << "{";
          for (const auto& cmd : cmds)
            os << sep << cmd, sep = "\n";
          return os << "}";
        }
      };
      const Piper& execute(const Commands& cmds) const {
        for (const auto& cmd : cmds)
          fputs((cmd + "\n").c_str(), file_.get());
        return *this;
      }

    private:
      std::unique_ptr<FILE, decltype(&pclose)> file_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
