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

#ifndef CepGen_Core_SteeredObject_h
#define CepGen_Core_SteeredObject_h

#include <algorithm>

#include "CepGen/Core/Steerable.h"

#define REGISTER_STEERED_OBJECT_CONTENT_TYPE   \
  __TYPE_ENUM(bool, map_bools_)                \
  __TYPE_ENUM(int, map_ints_)                  \
  __TYPE_ENUM(unsigned long long, map_ulongs_) \
  __TYPE_ENUM(double, map_dbls_)               \
  __TYPE_ENUM(std::string, map_strs_)          \
  __TYPE_ENUM(Limits, map_lims_)               \
  __TYPE_ENUM(ParametersList, map_params_)     \
  __TYPE_ENUM(std::vector<int>, map_vec_ints_) \
  __TYPE_ENUM(std::vector<Limits>, map_vec_lims_)

namespace cepgen {
  /// Base user-steerable object
  /// \tparam T Object steered (using its T::description() static member)
  template <typename T>
  class SteeredObject : public Steerable {
  public:
    SteeredObject() : Steerable(T::description().parameters()) {}  ///< Build a module
    explicit SteeredObject(const ParametersList& params) : Steerable(T::description().validate(params)) {}
    ~SteeredObject() override = default;

    /// Equality operator
    inline bool operator==(const SteeredObject& oth) const { return parameters() == oth.parameters(); }
    /// Inequality operator
    inline bool operator!=(const SteeredObject& oth) const { return !operator==(oth); }

    /// Module user-defined parameters
    inline const ParametersList& parameters() const override {
#define __TYPE_ENUM(type, map_name) \
  for (const auto& kv : map_name)   \
    params_.set(kv.first, kv.second);
      REGISTER_STEERED_OBJECT_CONTENT_TYPE
#undef __TYPE_ENUM
      return Steerable::parameters();
    }
    inline void setParameters(const ParametersList& params) override {
      if (params.empty())
        return;
      Steerable::setParameters(params);
#define __TYPE_ENUM(type, map_name) \
  for (const auto& kv : map_name)   \
    params_.fill(kv.first, kv.second);
      REGISTER_STEERED_OBJECT_CONTENT_TYPE
#undef __TYPE_ENUM
    }
    /// Set (documented) module parameters
    inline void setDescribedParameters(const ParametersList& params_orig) {
      const auto obj_keys = T::description().parameters().keys();
      if (obj_keys.empty())
        return;
      auto params = params_orig;
      for (const auto& key : params.keys())
        if (std::find(obj_keys.begin(), obj_keys.end(), key) == obj_keys.end())
          params.erase(key);
      setParameters(params);
    }

#define __TYPE_ENUM(type, map_name)                              \
  inline SteeredObject& add(const std::string& key, type& var) { \
    map_name.insert({key, var});                                 \
    map_name.at(key) = params_.operator[]<type>(key);            \
    return *this;                                                \
  }
    REGISTER_STEERED_OBJECT_CONTENT_TYPE
#undef __TYPE_ENUM

  private:
#define __TYPE_ENUM(type, map_name) std::unordered_map<std::string, type&> map_name;
    REGISTER_STEERED_OBJECT_CONTENT_TYPE
#undef __TYPE_ENUM
  };
}  // namespace cepgen

#undef REGISTER_STEERED_OBJECT_CONTENT_TYPE

#endif
