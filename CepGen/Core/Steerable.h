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

#ifndef CepGen_Core_Steerable_h
#define CepGen_Core_Steerable_h

#include "CepGen/Core/ParametersList.h"

namespace cepgen {
  /// Base runtime module object
  class Steerable {
  public:
    /// Build a module from its steering parameters
    explicit Steerable(const ParametersList& params) : params_(params) {}
    virtual ~Steerable() = default;

    /// Set all module parameters
    virtual void setParameters(const ParametersList& params) { params_ += params; }
    /// Module parameters
    virtual const ParametersList& parameters() const { return params_; }

    /// Retrieve a parameters as previously steered
    template <typename T>
    T steer(const std::string& key) const {
      return params_.get<T>(key);
    }
    /// Retrieve a recasted parameters as previously steered
    template <typename T, typename U>
    U steerAs(const std::string& key) const {
      return params_.getAs<T, U>(key);
    }

  protected:
    /// Module parameters
    mutable ParametersList params_;
  };
}  // namespace cepgen

#endif
