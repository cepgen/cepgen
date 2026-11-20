/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2025  Laurent Forthomme
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

#ifndef CepGen_Core_ParametersList_h
#define CepGen_Core_ParametersList_h

#include <map>
#include <string>
#include <unordered_map>

#include "CepGen/Utils/Limits.h"

/// Looper over the list of parameter containers handled by the ParametersList object
/// \note This can be edited to add an extra handled collection to this steering utility
///   e.g. add __TYPE_ENUM(typename,      // any C/C++ type name
///                        map_variable,  // ParametersList object attribute
///                        "human-readable name of parameter")
///        to the following list.
#define REGISTER_CONTENT_TYPE                                            \
  __TYPE_ENUM(bool, bool_values_, "bool")                                \
  __TYPE_ENUM(int, int_values_, "int")                                   \
  __TYPE_ENUM(unsigned long long, ulong_values_, "ulong")                \
  __TYPE_ENUM(double, dbl_values_, "float")                              \
  __TYPE_ENUM(std::string, str_values_, "str")                           \
  __TYPE_ENUM(Limits, lim_values_, "Limits")                             \
  __TYPE_ENUM(ParametersList, param_values_, "Params")                   \
  __TYPE_ENUM(std::vector<int>, vec_int_values_, "vint")                 \
  __TYPE_ENUM(std::vector<double>, vec_dbl_values_, "vfloat")            \
  __TYPE_ENUM(std::vector<std::string>, vec_str_values_, "vstr")         \
  __TYPE_ENUM(std::vector<Limits>, vec_lim_values_, "VLimits")           \
  __TYPE_ENUM(std::vector<ParametersList>, vec_param_values_, "VParams") \
  __TYPE_ENUM(std::vector<std::vector<double> >, vec_vec_dbl_values_, "vvfloat")

namespace cepgen {
  const auto MODULE_NAME = "mod_name";  ///< Indexing key for the module name
  /// Parameters container
  class ParametersList {
    /// Retrieve the default argument for a given variable type
    template <typename T>
    struct default_arg {
      static T get() { return T(); }  ///< Default variable argument
    };

  public:
    ParametersList() = default;
    ParametersList(const ParametersList&);  ///< Copy constructor
    ~ParametersList() = default;  // required for unique_ptr initialisation! avoids cleaning all individual objects
    ParametersList& operator=(const ParametersList&) = default;  ///< Assignment operator

    bool hasName() const;                                           ///< Does the parameters' list have a name key?
    std::string name(const std::string& default_value = "") const;  ///< Retrieve the module name if any
    ParametersList& setName(const std::string&);                    ///< Set the module name

    bool operator==(const ParametersList&) const;                                  ///< Equality operator
    bool operator!=(const ParametersList& oth) const { return !operator==(oth); }  ///< Inequality operator
    ParametersList diff(const ParametersList&) const;  ///< Diff with another parameters list ('mine' + 'theirs' keys)

    ParametersList& feed(const std::string&);  ///< Feed a control string to the list of parameters

    /// Check if a given parameter is handled in this list
    /// \param key Unique key for parameter
    template <typename T>
    bool has(const std::string& key) const;
    /// Erase a parameter with key
    /// \return Number of key-indexed values erased
    size_t erase(const std::string&);
    /// Erase a typed parameter with the key
    /// \return Number of key-indexed values erased
    template <typename T>
    size_t erase(const std::string&);

    /// Fill a variable with the key content if exists
    /// \param[in] key Unique key for the parameter
    /// \param[out] value Object/variable to fill with parameter value
    /// \note If the key is stored with another type, its value is losslessly converted (see cast)
    /// \throw Exception if the key is stored with a type that cannot be converted into the requested one
    template <typename T>
    const ParametersList& fill(const std::string& key, T& value) const {
      if (has<T>(key) || !typeOf(key).empty())
        value = get<T>(key);
      return *this;
    }
    /// Get a parameter value
    /// \param[in] key Unique key for the parameter
    /// \param[in] def Default parameters values if parameter is not contained
    /// \note If the key is stored with another type, its value is losslessly converted (see cast)
    /// \throw Exception if the key is stored with a type that cannot be converted into the requested one
    template <typename T>
    T get(const std::string& key, const T& def = default_arg<T>::get()) const;
    /// Fill a variable with a parameter value, converted into the requested type if this conversion is lossless
    /// \param[in] key Unique key for the parameter
    /// \param[out] value Object/variable to fill with parameter value, left untouched if the conversion failed
    /// \return False if the key is not found, or stored with a type not convertible into the requested one
    template <typename T>
    bool cast(const std::string& key, T& value) const {
      if (has<T>(key)) {
        value = get<T>(key);
        return true;
      }
      return convert<T>(key, value);
    }
    /// Human-readable name of the type a parameter is stored with (e.g. "float", "Limits")
    /// \return Empty string if the key is not found
    std::string typeOf(const std::string& key) const;
    /// Get a recast parameter value
    /// \tparam T Base parameter type
    /// \tparam U Type to recast the parameter into
    /// \param[in] key Unique key for the parameter
    /// \param[in] def Default parameters values if parameter is not contained
    template <typename T, typename U>
    U getAs(const std::string& key, const U& def = default_arg<U>::get()) const {
      return static_cast<U>(get<T>(key, static_cast<T>(def)));
    }
    /// Reference to a parameter value
    /// \param[in] key Unique key for parameter
    template <typename T>
    T& operator[](const std::string& key);
    /// Set a parameter value
    /// \note Any value stored with the same key but another type is erased: a key is only held with one type
    template <typename T>
    ParametersList& set(const std::string&, const T&);
    /// Set a recast parameter value
    /// \tparam T Base parameter type
    /// \tparam U Type to recast the parameter into
    /// \param[in] key Unique key for the parameter
    /// \param[in] value Value to set the parameter
    template <typename T, typename U>
    ParametersList& setAs(const std::string& key, const U& value) {
      return set<T>(key, static_cast<T>(value));
    }
    ParametersList& rename(const std::string&, const std::string&);  ///< Rename the key to a parameter value
    ParametersList& operator+=(const ParametersList& oth);           ///< Concatenate two parameters containers
    ParametersList operator+(const ParametersList& other) const;     ///< Concatenation of two-parameter containers
    bool empty() const;                                              ///< Is the list empty?

    template <typename T>
    std::vector<std::string> keysOf() const;  ///< List of keys for one type in this list of parameters
    /// List of keys handled in this list of parameters
    /// \param[in] name_key Include the name variable?
    std::vector<std::string> keys(bool name_key = true) const;
    /// Get a string-converted version of a value
    /// \param[in] key Unique key for the parameter
    /// \param[in] wrap Encapsulate the value with type()
    std::string getString(const std::string& key, bool wrap = false) const;
    /// Get a string-converted version of the module name if any
    /// \param[in] wrap Encapsulate the value with type()
    inline std::string getNameString(bool wrap = false) const { return getString(MODULE_NAME, wrap); }
    std::string serialise() const;  ///< Serialise a parameter collection into a parseable string

    friend std::ostream& operator<<(std::ostream&,
                                    const ParametersList&);  ///< Human-readable version of a parameter container
    const ParametersList& print(std::ostream&) const;        ///< Debugging-like printout of a parameter container
    std::string print(bool compact = false) const;           ///< Normal printout of a parameter container

  private:
    /// Lossless conversion of a value stored with another type than the requested one
    template <typename T>
    bool convert(const std::string&, T&) const {
      return false;
    }
    /// Is the key holding an empty vector (of any type)? An empty vector can be converted into any vector type
    bool isEmptyVector(const std::string&) const;
    /// Erase all values with a given key but stored with another type than the requested one
    template <typename T>
    void eraseOtherTypes(const std::string&);

    std::map<std::string, ParametersList> param_values_;
    std::unordered_map<std::string, bool> bool_values_;
    std::unordered_map<std::string, int> int_values_;
    std::unordered_map<std::string, unsigned long long> ulong_values_;
    std::unordered_map<std::string, double> dbl_values_;
    std::unordered_map<std::string, std::string> str_values_;
    std::unordered_map<std::string, Limits> lim_values_;
    std::unordered_map<std::string, std::vector<int> > vec_int_values_;
    std::unordered_map<std::string, std::vector<double> > vec_dbl_values_;
    std::unordered_map<std::string, std::vector<std::string> > vec_str_values_;
    std::unordered_map<std::string, std::vector<Limits> > vec_lim_values_;
    std::unordered_map<std::string, std::vector<ParametersList> > vec_param_values_;
    std::unordered_map<std::string, std::vector<std::vector<double> > > vec_vec_dbl_values_;
  };

/// Implement all setters and getters for a given type
#define __TYPE_ENUM(type, map, name)                                          \
  template <>                                                                 \
  bool ParametersList::has<type>(const std::string&) const;                   \
  template <>                                                                 \
  type ParametersList::get<type>(const std::string&, const type&) const;      \
  template <>                                                                 \
  type& ParametersList::operator[]<type>(const std::string&);                 \
  template <>                                                                 \
  ParametersList& ParametersList::set<type>(const std::string&, const type&); \
  template <>                                                                 \
  std::vector<std::string> ParametersList::keysOf<type>() const;              \
  template <>                                                                 \
  size_t ParametersList::erase<type>(const std::string&);                     \
  template <>                                                                 \
  bool ParametersList::convert<type>(const std::string&, type&) const;
  REGISTER_CONTENT_TYPE
#undef __TYPE_ENUM
}  // namespace cepgen

#endif
