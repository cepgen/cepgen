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

#ifndef CepGen_Physics_Utils_h
#define CepGen_Physics_Utils_h

namespace cepgen {
  namespace utils {
    /// Compute the diffractive mass from virtuality/Bjorken x
    double mX2(double xbj, double q2, double mp2);
    /// Compute Bjorken x from virtuality/diffractive mass
    double xBj(double q2, double mp2, double mx2);
    /// Compute the virtuality from Bjorken x/diffractive mass
    double q2(double xbj, double mp2, double mx2);
    /// Compute energy from mass and emitted mass
    double energyFromW(double w, double mp2, double m2);
    namespace kt {
      /// Compute the diffractive mass from longitudinal loss/transverse virtuality/virtuality
      double mX2(double x, double kt2, double q2, double mi2);
      /// Compute the virtuality from longitudinal loss/transverse virtuality/diffractive mass
      double q2(double x, double kt2, double mi2, double mf2 = -1.);
    }  // namespace kt
  }    // namespace utils
}  // namespace cepgen

#endif
