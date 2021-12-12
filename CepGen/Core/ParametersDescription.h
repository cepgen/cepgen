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

#ifndef CepGen_Core_ParametersDescription_h
#define CepGen_Core_ParametersDescription_h

#include "CepGen/Core/ParametersList.h"

namespace cepgen {
  /// A description object for parameters collection
  class ParametersDescription : private ParametersList {
  public:
    /// Build the description of a parameters collection object
    /// \param[in] mod_name Module name (where applicable)
    explicit ParametersDescription(const std::string& mod_name = "");
    /// Build the (empty) description of a parameters collection object from its definition
    explicit ParametersDescription(const ParametersList& params);
    /// Copy constructor
    ParametersDescription(const ParametersDescription&) = default;

    /// Does a description of this parameter (or parameters collection) exist?
    bool empty() const;
    /// Concatenate another description to this one
    ParametersDescription& operator+=(const ParametersDescription&);
    /// Human-readable description
    friend std::ostream& operator<<(std::ostream&, const ParametersDescription&);
    /// Set the module name for this parameter (or parameters collection)
    template <typename T>
    ParametersDescription& setName(const T& name) {
      add<T>(ParametersList::MODULE_NAME, name);
      return *this;
    }
    /// Set the description of this parameter (or parameters collection)
    ParametersDescription& setDescription(const std::string& descr);
    /// Description of this parameter (or parameters collection)
    const std::string& description() const { return mod_descr_; }
    /// Add the description to a new parameter
    template <typename T>
    ParametersDescription& add(const std::string& name, const T& def) {
      obj_descr_[name] = ParametersDescription();
      ParametersList::set<T>(name, def);
      return obj_descr_[name];
    }
    /// Add the description to a collection of ParametersList objects
    ParametersDescription& addParametersDescriptionVector(const std::string&, const ParametersDescription&);
    /// Human-readable description of all parameters and their default value
    std::string describe(size_t offset = 0) const;
    /// List of parameters associated to this description object
    ParametersList& parameters();
    /// List of parameters associated to this description object
    const ParametersList& parameters() const;
    /// Get the description of a sub-object
    const ParametersDescription& get(const std::string&) const;

    /// Parameter type
    enum struct Type { Value, Parameters, Module };
    /// Get the type of parameter considered
    Type type() const;
    /// Validate a set of used-steered parameters
    void validate(const ParametersList&) const;

  private:
    std::string mod_descr_;
    std::map<std::string, ParametersDescription> obj_descr_;
  };
  /// Add the description to a new sub-description (aka ParametersList) object
  template <>
  ParametersDescription& ParametersDescription::add(const std::string&, const ParametersDescription&);
  /// Disable the addition of a ParametersList object to this description
  template <>
  ParametersDescription& ParametersDescription::add(const std::string&, const ParametersList&);
}  // namespace cepgen

#endif
