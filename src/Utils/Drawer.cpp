/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#include <bitset>

#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Message.h"

using namespace cepgen::utils;
using namespace std::string_literals;

Drawer::Drawer(const ParametersList& params) : NamedModule(params) {}

Drawer::Mode operator|(const Drawer::Mode& lhs, const Drawer::Mode::value_t& rhs) {
  std::bitset<16> mod1(static_cast<int>(lhs.value())), mod2(static_cast<int>(rhs));
  return Drawer::Mode(static_cast<Drawer::Mode::value_t>((mod1 | mod2).to_ulong()));
}

namespace cepgen::utils {
  bool operator&(const Drawer::Mode& lhs, const Drawer::Mode::value_t& rhs) {
    return static_cast<int>(lhs.value_) & static_cast<int>(rhs);
  }

  std::ostream& operator<<(std::ostream& os, const Drawer::Mode& mode) {
    if (mode.value() == Drawer::Mode::none)
      return os << "none";
    std::string sep;
    if (mode & Drawer::Mode::logx)
      os << sep << "logx"s, sep = "|";
    if (mode & Drawer::Mode::logy)
      os << sep << "logy"s, sep = "|";
    if (mode & Drawer::Mode::logz)
      os << sep << "logz"s, sep = "|";
    if (mode & Drawer::Mode::nostack)
      os << sep << "nostack"s, sep = "|";
    if (mode & Drawer::Mode::grid)
      os << sep << "grid", sep = "|";
    if (mode & Drawer::Mode::col)
      os << sep << "col", sep = "|";
    if (mode & Drawer::Mode::cont)
      os << sep << "cont", sep = "|";
    if (mode & Drawer::Mode::ratio)
      os << sep << "ratio";
    return os;
  }
}  // namespace cepgen::utils

namespace cepgen {
  Drawer::Mode operator|(const Drawer::Mode::value_t& lhs, const Drawer::Mode::value_t& rhs) {
    const std::bitset<16> mod1(static_cast<int>(lhs)), mod2(static_cast<int>(rhs));
    return Drawer::Mode((mod1 | mod2).to_ulong());
  }
}  // namespace cepgen

Drawer::Mode& operator|=(Drawer::Mode& one, const Drawer::Mode::value_t& oth) {
  one = one | oth;
  return one;
}
