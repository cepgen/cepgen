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

#ifndef CepGen_Version_h
#define CepGen_Version_h

#include <string>

namespace cepgen {
  /// Collection of CepGen version information handlers
  struct version {
    /// CepGen version
    static const std::string tag;
    /// CepGen detailed version
    static const std::string extended;
    /// CepGen banner
    static const std::string banner;
  };
}  // namespace cepgen

#endif
