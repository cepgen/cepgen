/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#ifndef CepGen_Utils_Caller_h
#define CepGen_Utils_Caller_h

#include <iosfwd>
#include <sstream>
#include <string>
#include <vector>

namespace cepgen::utils {
  /// External command piping utility
  class Caller {
  public:
    Caller();
    ~Caller();
    /// Start a logged call command
    /// \param[in] commands Command path for the session
    static std::string call(const std::vector<std::string>& commands);
    /// Start a logged call command
    /// \param[in] command Command path for the session
    static std::string call(const std::string& command);

    std::string output() const;  ///< Retrieve the (potential) output from the command
    std::string error() const;   ///< Retrieve the (potential) error stream from the command

  private:
    std::stringstream os_cout_, os_cerr_;
    std::streambuf *oldcout_{nullptr}, *oldcerr_{nullptr};
  };
}  // namespace cepgen::utils

#endif
