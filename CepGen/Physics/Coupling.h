/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

#ifndef CepGen_Physics_Coupling_h
#define CepGen_Physics_Coupling_h

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  /// A generic \f$\alpha_S(Q^2)\f$ evaluation algorithm
  class Coupling : public NamedModule<Coupling, std::string> {
  public:
    /// Build an \f$\alpha_{S,EM}\f$ interpolator object
    explicit Coupling(const ParametersList& params) : NamedModule(params) {}
    virtual ~Coupling() {}
    virtual double operator()(double q) const = 0;  ///< Compute \f$\alpha_{S,EM}\f$ for a given \f$Q\f$
  };
}  // namespace cepgen

#endif
