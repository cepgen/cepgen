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

#ifndef CepGen_Utils_Physics_h
#define CepGen_Utils_Physics_h

namespace cepgen {
  namespace utils {
    /// Compute the diffractive mass from virtuality/Bjorken x
    double mX2(double xbj, double q2, double mp2);
    /// Compute Bjorken x from virtuality/diffractive mass
    double xBj(double q2, double mp2, double mx2);
  }  // namespace utils
}  // namespace cepgen

#endif
