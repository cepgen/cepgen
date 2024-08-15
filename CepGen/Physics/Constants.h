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

#ifndef CepGen_Physics_Constants_h
#define CepGen_Physics_Constants_h

#include <cmath>

/// List of physical constants useful that may be used for the matrix element definition
namespace cepgen::constants {
  /// Electromagnetic coupling constant \f$\alpha_{\rm em}=\frac{e^2}{4\pi\epsilon_0\hbar c}\f$
  constexpr double ALPHA_EM = 1. / 137.035999;
#if !defined(__CINT__) && !defined(__CLING__)
  constexpr double G_EM_SQ = 4. * M_PI * ALPHA_EM;  ///< Electromagnetic charge (~0.303 in natural units)
#endif
  constexpr double ALPHA_QCD = 0.1184;  ///< Strong coupling constant \f$\alpha_{\rm QCD} at the Z pole\f$
  /// Conversion factor between GeV\f$^{-2}\f$ and barn
  /// i.e. \f$\hbar^2 c^2\f$ in GeV\f$^{-2}\f$.
  constexpr double GEVM2_TO_PB = 0.389351824e9;
  /// Conversion factor between GeV\f$^{-1}\f$ and metres
  /// i.e. \f$\hbar c\f$ in GeV\f$^{-1}\f$.
  constexpr double GEVM1_TO_M = 0.1973269804e-15;
  constexpr double G_F = 1.1663787e-5;  ///< Fermi coupling constant, in GeV\f$^{-2}\f$
}  // namespace cepgen::constants

#endif
