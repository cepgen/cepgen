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

#ifndef CepGen_Core_ParametersDescription_h
#define CepGen_Core_ParametersDescription_h

#include "CepGen/Core/ParametersList.h"

namespace cepgen {
  /// A description object for parameters collection
  class ParametersDescription : private ParametersList {
  public:
    /// Build the description of a parameters collection object
    /// \param[in] mod_key Module name (where applicable)
    explicit ParametersDescription(const std::string& mod_key = "");
    /// Build the (empty) description of a parameters collection object from its definition
    explicit ParametersDescription(const ParametersList& params);

    /// Does a description of this parameter (or parameters collection) exist?
    bool empty() const;
    /// Concatenate another description to this one
    ParametersDescription& operator+=(const ParametersDescription&);
    /// Human-readable description
    friend std::ostream& operator<<(std::ostream&, const ParametersDescription&);
    /// Set the module name for this parameter (or parameters collection)
    template <typename I>
    ParametersDescription& setKey(const I& key) {
      mod_key_ = std::to_string(key);
      return *this;
    }
    /// Module name for this parameter
    const std::string& key() const { return mod_key_; }
    /// Set the description of this parameter (or parameters collection)
    ParametersDescription& setDescription(const std::string& descr);
    /// Description of this parameter (or parameters collection)
    const std::string& description() const { return mod_descr_; }
    /// This parameters is a collection of sub-parameters
    ParametersDescription& setParametersVector(bool pv = true) {
      is_vec_params_ = pv;
      return *this;
    }
    /// Add the description to a new parameter
    template <typename T>
    ParametersDescription& add(const std::string& name, const T& def) {
      if (obj_descr_.count(name) == 0)
        // only add a new, empty description if not yet described
        // (allows to ensure previous descriptions are not discarded)
        obj_descr_[name] = ParametersDescription();
      ParametersList::set<T>(name, def);
      return obj_descr_[name];
    }
    /// Add a recast definition to a new parameter
    template <typename T, typename U>
    inline ParametersDescription& addAs(const std::string& name, const U& def) {
      return add<T>(name, static_cast<T>(def));
    }
    /// Set the module name
    inline ParametersDescription& setName(const std::string& name) {
      ParametersList::setName(name);
      return *this;
    }
    /// Add the description to a collection of ParametersList objects
    ParametersDescription& addParametersDescriptionVector(const std::string&,
                                                          const ParametersDescription&,
                                                          const std::vector<ParametersList>& def = {});
    std::string describe(size_t offset = 0) const;  ///< Human-readable description of parameters and their default value
    ParametersList& parameters();                   ///< List of parameters associated to this description object
    const ParametersList& parameters() const;       ///< List of parameters associated to this description object
    bool has(const std::string&) const;             ///< Ensure the description exists
    const ParametersDescription& get(const std::string&) const;  ///< Get the description of a sub-object
    ParametersList validate(const ParametersList&) const;        ///< Validate a set of used-steered parameters
    ParametersDescription steer(const ParametersList&) const;    ///< Set parameters value for this description object

    enum struct Type { Value, Parameters, Module, ParametersVector };  ///< Parameter type
    friend std::ostream& operator<<(std::ostream&, const Type&);       ///< Human-readable description of parameter type
    Type type() const;                                                 ///< Get the type of parameter considered

    /// A collection of valid values for a given parameter
    class ParameterValues {
    public:
      ParameterValues() {}

      /// Short printout of allowed parameter values
      friend std::ostream& operator<<(std::ostream&, const ParameterValues&);

      bool empty() const;  ///< Check if a parameter has a limited set of allowed values

      ParameterValues& allow(int, const std::string& = "");                 ///< Allow an integer value for a parameter
      ParameterValues& allow(const std::string&, const std::string& = "");  ///< Allow a string value for a parameter
      std::map<std::string, std::string> allowed() const;  ///< Helper list of allowed values (all types) for a parameter

      bool validate(int) const;                 ///< Check if an integer value is allowed for a parameter
      bool validate(const std::string&) const;  ///< Check if a string value is allowed for a parameter

    private:
      std::map<int, std::string> int_vals_;
      std::map<std::string, std::string> str_vals_;
    };

    inline const ParameterValues& values() const { return obj_values_; }  ///< Possible values for a parameter
    inline ParameterValues& values() { return obj_values_; }              ///< Possible values for a parameter

  private:
    std::string mod_key_, mod_descr_;
    bool is_vec_params_{false};
    std::map<std::string, ParametersDescription> obj_descr_;
    ParameterValues obj_values_;
  };
  template <>
  ParametersDescription& ParametersDescription::setKey<std::string>(const std::string&);
  /// Add the description to a new sub-description (aka ParametersList) object
  template <>
  ParametersDescription& ParametersDescription::add(const std::string&, const ParametersDescription&);
  /// Disable the addition of a ParametersList object to this description
  template <>
  ParametersDescription& ParametersDescription::add(const std::string&, const ParametersList&);
}  // namespace cepgen

#endif
