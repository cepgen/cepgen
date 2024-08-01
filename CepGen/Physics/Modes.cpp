/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

#include "CepGen/Physics/Modes.h"

namespace cepgen::mode {
  std::ostream& operator<<(std::ostream& os, const Kinematics& pm) {
    switch (pm) {
      case Kinematics::invalid:
        return os << "{invalid}";
      case Kinematics::ElasticElastic:
        return os << "elastic/elastic";
      case Kinematics::InelasticElastic:
        return os << "inelastic/elastic";
      case Kinematics::ElasticInelastic:
        return os << "elastic/inelastic";
      case Kinematics::InelasticInelastic:
        return os << "inelastic/inelastic";
    }
    return os;
  }
}  // namespace cepgen::mode
