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

#ifndef CepGen_Integration_Integrand_h
#define CepGen_Integration_Integrand_h

#include <vector>

namespace cepgen {
  /// An integrand wrapper placeholder
  class Integrand {
  public:
    Integrand() = default;
    virtual ~Integrand() = default;

    virtual double eval(const std::vector<double>&) = 0;  ///< Compute the integrand for a given coordinates set
    virtual size_t size() const = 0;                      ///< Phase space dimension
    virtual bool hasProcess() const { return false; }     ///< Does this integrand also contain a process object?
  };
}  // namespace cepgen

#endif
