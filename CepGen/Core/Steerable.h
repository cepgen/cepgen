/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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

#include "CepGen/Core/ParametersDescription.h"

namespace cepgen {
  /// Base runtime module object
  class Steerable {
  public:
    explicit Steerable(const ParametersList&);  ///< Build a module from its steering parameters
    virtual ~Steerable() = default;

    inline virtual void setParameters(const ParametersList& params) { params_ += params; }  ///< Set module parameters
    inline virtual const ParametersList& parameters() const { return params_; }             ///< Module parameters

    /// Description of all object parameters
    static inline ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.setDescription("Virtual, base steerable object");
      return desc;
    }

  protected:
    /// Retrieve a parameters as previously steered
    template <typename T>
    inline T steer(const std::string& key) const {
      return params_.get<T>(key);
    }
    /// Retrieve a recasted parameters as previously steered
    template <typename T, typename U>
    inline U steerAs(const std::string& key) const {
      return params_.getAs<T, U>(key);
    }
    /// Retrieve module name from parameters
    template <typename T>
    inline T steerName() const {
      return steer<T>(MODULE_NAME);
    }
    std::string steerPath(const std::string& key) const;  ///< Retrieve a path from common search paths

    mutable ParametersList params_;  ///< Module parameters
  };
}  // namespace cepgen

#endif
