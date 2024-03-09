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

#include <functional>  // for std::reference_wrapper

#include "CepGen/Core/Steerable.h"

#define REGISTER_TYPE(type, coll)                                      \
public:                                                                \
  inline SteeredObject& add(const std::string& key, type& var) {       \
    coll.insert({key, std::ref(var)});                                 \
    coll.at(key).get() = params_.operator[]<type>(key);                \
    return *this;                                                      \
  }                                                                    \
                                                                       \
private:                                                               \
  std::unordered_map<std::string, std::reference_wrapper<type> > coll; \
  static_assert(true, "")

namespace cepgen {
  /// Base user-steerable object
  /// \tparam T the type of object to be steered (using its T::description() static member)
  template <typename T>
  class SteeredObject : public Steerable {
  public:
    /// Build a module
    inline SteeredObject() : Steerable(T::description().parameters()) {}
    explicit inline SteeredObject(const ParametersList& params) : Steerable(T::description().validate(params)) {}
    virtual ~SteeredObject() = default;

    /// Equality operator
    inline bool operator==(const SteeredObject& oth) const { return parameters() == oth.parameters(); }
    /// Inequality operator
    inline bool operator!=(const SteeredObject& oth) const { return !operator==(oth); }

    REGISTER_TYPE(bool, map_bools_);
    REGISTER_TYPE(int, map_ints_);
    REGISTER_TYPE(unsigned long long, map_ulongs_);
    REGISTER_TYPE(double, map_dbls_);
    REGISTER_TYPE(std::string, map_strs_);
    REGISTER_TYPE(Limits, map_lims_);
    REGISTER_TYPE(ParametersList, map_params_);
    REGISTER_TYPE(std::vector<int>, map_vints_);

  public:
    /// Module user-defined parameters
    inline const ParametersList& parameters() const override {
      const auto set = [this](auto& map) {
        for (const auto& kv : map)
          params_.set(kv.first, kv.second.get());
      };
      set(map_bools_);
      set(map_ints_);
      set(map_ulongs_);
      set(map_dbls_);
      set(map_strs_);
      set(map_lims_);
      set(map_params_);
      set(map_vints_);
      return params_;
    }
    virtual inline void setParameters(const ParametersList& params) override {
      if (params.empty())
        return;
      Steerable::setParameters(params);
      const auto fill = [this](auto& map) {
        for (const auto& kv : map)
          params_.fill(kv.first, kv.second.get());
      };
      fill(map_bools_);
      fill(map_ints_);
      fill(map_ulongs_);
      fill(map_dbls_);
      fill(map_strs_);
      fill(map_lims_);
      fill(map_params_);
      fill(map_vints_);
    }
  };
}  // namespace cepgen

#undef REGISTER_TYPE

#endif
