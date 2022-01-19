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

#include <gsl/gsl_errno.h>

#include <cmath>
#include <iomanip>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Drawable.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    Drawable::Drawable(const std::string& name, const std::string& title) : name_(name), title_(title) {}

    Drawable::Drawable(const Drawable& oth) : xlabel_(oth.xlabel_), ylabel_(oth.ylabel_) {}

    Drawer::Mode operator|(const Drawer::Mode& lhs, const Drawer::Mode& rhs) {
      std::bitset<7> mod1((int)lhs), mod2((int)rhs);
      return (Drawer::Mode)(mod1 | mod2).to_ulong();
    }

    bool operator&(const Drawer::Mode& lhs, const Drawer::Mode& rhs) {
      //return ((int)lhs > (int)rhs) - ((int)lhs < (int)rhs);
      return (int)lhs & (int)rhs;
    }
  }  // namespace utils
}  // namespace cepgen
