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

#include <bitset>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"

namespace cepgen {
  namespace utils {
    Drawer::Drawer(const ParametersList& params) : NamedModule(params) {}

    Drawer::Mode operator|(const Drawer::Mode& lhs, const Drawer::Mode& rhs) {
      std::bitset<16> mod1((int)lhs.value_), mod2((int)rhs.value_);
      Drawer::Mode out;
      out.value_ = (Drawer::Mode::value_t)(mod1 | mod2).to_ulong();
      return out;
    }

    Drawer::Mode& operator|=(Drawer::Mode& one, const Drawer::Mode& oth) {
      one = one | oth.value_;
      return one;
    }

    bool operator&(const Drawer::Mode& lhs, const Drawer::Mode::value_t& rhs) {
      //return ((int)lhs > (int)rhs) - ((int)lhs < (int)rhs);
      return (int)lhs.value_ & (int)rhs;
    }

    std::ostream& operator<<(std::ostream& os, const Drawer::Mode& mode) {
      if (mode.value_ == Drawer::Mode::none)
        return os << "none";
      std::string sep;
      if (mode & Drawer::Mode::logx)
        os << sep << "logx", sep = "|";
      if (mode & Drawer::Mode::logy)
        os << sep << "logy", sep = "|";
      if (mode & Drawer::Mode::logz)
        os << sep << "logz", sep = "|";
      if (mode & Drawer::Mode::nostack)
        os << sep << "nostack", sep = "|";
      if (mode & Drawer::Mode::grid)
        os << sep << "grid", sep = "|";
      if (mode & Drawer::Mode::col)
        os << sep << "col", sep = "|";
      if (mode & Drawer::Mode::cont)
        os << sep << "cont", sep = "|";
      return os;
    }
  }  // namespace utils
}  // namespace cepgen
