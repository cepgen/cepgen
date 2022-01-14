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

#ifndef CepGen_Utils_Drawer_h
#define CepGen_Utils_Drawer_h

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  namespace utils {
    class Drawable;
    /// A generic drawing utilitary
    class Drawer : public NamedModule<std::string> {
    public:
      /// Build an \f$\alpha_{S,EM}\f$ interpolator object
      explicit Drawer(const ParametersList& params) : NamedModule(params) {}
      /// Compute \f$\alpha_{S,EM}\f$ for a given \f$Q\f$
      virtual const Drawer& draw(const Drawable&) const = 0;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
