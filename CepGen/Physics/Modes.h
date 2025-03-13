/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#ifndef CepGen_Physics_Modes_h
#define CepGen_Physics_Modes_h

#include <iosfwd>

/// Collection of enums for the definition of process mode
namespace cepgen::mode {
  /// Type of scattering
  enum class Kinematics {
    invalid = 0,
    ElasticElastic = 1,     ///< proton-proton elastic case
    ElasticInelastic = 2,   ///< proton-proton single-dissociative (or inelastic-elastic) case
    InelasticElastic = 3,   ///< proton-proton single-dissociative (or elastic-inelastic) case
    InelasticInelastic = 4  ///< proton-proton double-dissociative case
  };
  /// Human-readable format of a process mode (elastic/dissociative parts)
  std::ostream& operator<<(std::ostream&, const Kinematics&);
}  // namespace cepgen::mode

#endif
