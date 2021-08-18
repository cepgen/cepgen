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

#ifndef CepGen_Physics_MCDFileParser_h
#define CepGen_Physics_MCDFileParser_h

#include <string>
#include <unordered_map>

namespace pdg {
  /// A MCD files parsing module
  class MCDFileParser {
  public:
    MCDFileParser() = default;
    /// Parse an external MCD file and retrieve all particles definition
    static void parse(const std::string& path);

  private:
    static constexpr size_t PDG_BEG = 1, PDG_END = 33;
    static constexpr size_t MASS_BEG = 33, MASS_END = 70;
    static constexpr size_t WIDTH_BEG = 70, WIDTH_END = 107;
    static constexpr size_t AUX_BEG = 107;
    static const std::unordered_map<std::string, short> MAP_CHARGE_STR;
  };
}  // namespace pdg

#endif
