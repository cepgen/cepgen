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

#include "CepGen/Utils/Physics.h"

namespace cepgen {
  namespace utils {
    double mX2(double xbj, double q2, double mp2) { return mp2 + q2 * (1. - xbj) / xbj; }

    double xBj(double q2, double mp2, double mx2) { return q2 / (q2 - mp2 + mx2); }
  }  // namespace utils
}  // namespace cepgen
