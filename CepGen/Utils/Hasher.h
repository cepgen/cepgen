/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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

#ifndef CepGen_Utils_Hasher_h
#define CepGen_Utils_Hasher_h

#include <cstddef>
#include <functional>

namespace cepgen::utils {
  /// A hasher table for a given structure
  template <class T, bool>
  struct Hasher {
    /// Hash a generic table
    inline size_t operator()(const T& t) const { return std::hash<T>()(t); }
  };
  /// A hasher table for a given structure
  template <class T>
  struct Hasher<T, true> {
    /// Hash a structure-indexed table
    inline size_t operator()(const T& t) {
      using enumType = typename std::underlying_type<T>::type;
      return std::hash<enumType>()(static_cast<enumType>(t));
    }
  };
  /// A hasher table for an enumeration
  template <class T>
  struct EnumHash {
    /// Hash an enumerator-indexed table
    inline size_t operator()(const T& t) const { return Hasher<T, std::is_enum<T>::value>()(t); }
  };
}  // namespace cepgen::utils

#endif
