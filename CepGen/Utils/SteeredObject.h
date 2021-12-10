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

#ifndef CepGen_Utils_SteeredObject_h
#define CepGen_Utils_SteeredObject_h

#include "CepGen/Utils/Steerable.h"

#define REGISTER_TYPE(type, coll)                                \
public:                                                          \
  inline SteeredObject& add(const std::string& key, type& var) { \
    coll[key] = &var;                                            \
    var = params_.operator[]<type>(key);                         \
    return *this;                                                \
  }                                                              \
                                                                 \
private:                                                         \
  std::unordered_map<std::string, type*> coll;

namespace cepgen {
  /// Base user-steerable object
  template <typename T>
  class SteeredObject : public Steerable {
  public:
    /// Build a module
    SteeredObject() : Steerable(T::description().parameters()) {}
    SteeredObject(const ParametersList& params) : Steerable(T::description().parameters() + params) {}
    virtual ~SteeredObject() = default;

    /// Description of all object parameters
    static inline ParametersDescription description() {
      auto desc = Steerable::description();
      desc.setDescription("Virtual, base steerable object");
      return desc;
    }

    /// Module user-defined parameters
    inline const ParametersList& parameters() const override {
      for (const auto& kv : map_bools_)
        params_.set<bool>(kv.first, *kv.second);
      for (const auto& kv : map_ints_)
        params_.set<int>(kv.first, *kv.second);
      for (const auto& kv : map_dbls_)
        params_.set<double>(kv.first, *kv.second);
      for (const auto& kv : map_strs_)
        params_.set<std::string>(kv.first, *kv.second);
      return params_;
    }

    void setParameters(const ParametersList& params) override { params_ = T::description().parameters() + params; }

    REGISTER_TYPE(bool, map_bools_)
    REGISTER_TYPE(int, map_ints_)
    REGISTER_TYPE(double, map_dbls_)
    REGISTER_TYPE(std::string, map_strs_)
  };
}  // namespace cepgen

#undef REGISTER_TYPE

#endif
