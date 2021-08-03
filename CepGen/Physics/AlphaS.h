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

#ifndef CepGen_Physics_AlphaS_h
#define CepGen_Physics_AlphaS_h

#include "CepGen/Modules/ModuleFactory.h"
#include "CepGen/Modules/NamedModule.h"

#define REGISTER_ALPHAS_MODULE(name, obj)                                  \
  namespace cepgen {                                                       \
    struct BUILDERNM(obj) {                                                \
      BUILDERNM(obj)() { AlphaSFactory::get().registerModule<obj>(name); } \
    };                                                                     \
    static const BUILDERNM(obj) g##obj;                                    \
  }

namespace cepgen {
  /// A generic \f$\alpha_S(Q^2)\f$ evaluation algorithm
  class AlphaS : public NamedModule<> {
  public:
    /// Build an \f$\alpha_S\f$ interpolator object
    AlphaS(const ParametersList& params) : NamedModule(params) {}
    /// Compute \f$\alpha_S\f$ for a given \f$Q^2\f$
    virtual double operator()(double q) const = 0;
  };
  /// An alpha(S) evolution algorithms factory
  typedef ModuleFactory<AlphaS> AlphaSFactory;
}  // namespace cepgen

#endif
