/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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

/// Looper over the list of parameters containers handled by the ParametersList object
/// \note This can be edited to add an extra handled collection to this steering utility
///   e.g. add __TYPE_ENUM(typename,      // any C/C++ type name
///                        map_variable,  // ParametersList object private attribute
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
  const char* const MODULE_NAME = "mod_name";  ///< Indexing key for the module name
  /// Parameters container
  class ParametersList {
  private:
    /// Retrieve the default argument for a given variable type
    template <typename T>
    struct default_arg {
      static inline T get() { return T(); }  ///< Default variable argument
    };

  public:
    ParametersList() = default;
    ParametersList(const ParametersList&);  ///< Copy constructor
    ~ParametersList() {}  // required for unique_ptr initialisation! avoids cleaning all individual objects
    ParametersList& operator=(const ParametersList&) = default;  ///< Assignment operator

    bool hasName() const;                                 ///< Does the parameters list have a name key?
    std::string name(const std::string& def = "") const;  ///< Retrieve the module name if any
    ParametersList& setName(const std::string&);          ///< Set the module name

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
    /// Erase a typed parameter with key
    /// \return Number of key-indexed values erased
    template <typename T>
    size_t erase(const std::string&);

    /// Fill a variable with the key content if exists
    /// \param[in] key Unique key for parameter
    /// \param[out] value Object/variable to be filled with parameter value
    template <typename T>
    inline const ParametersList& fill(const std::string& key, T& value) const {
      if (has<T>(key))
        value = get<T>(key);
      return *this;
    }
    /// Get a parameter value
    /// \param[in] key Unique key for parameter
    /// \param[in] def Default parameters value if parameter is not contained
    template <typename T>
    T get(const std::string& key, const T& def = default_arg<T>::get()) const;
    /// Get a recast parameter value
    /// \tparam T Base type of the parameter
    /// \tparam U Type to recast the parameter into
    /// \param[in] key Unique key for parameter
    /// \param[in] def Default parameters value if parameter is not contained
    template <typename T, typename U>
    inline U getAs(const std::string& key, const U& def = default_arg<U>::get()) const {
      return static_cast<U>(get<T>(key, static_cast<T>(def)));
    }
    /// Reference to a parameter value
    /// \param[in] key Unique key for parameter
    template <typename T>
    T& operator[](const std::string& key);
    template <typename T>
    ParametersList& set(const std::string&, const T&);  ///< Set a parameter value
    /// Set a recast parameter value
    /// \tparam T Base type of the parameter
    /// \tparam U Type to recast the parameter into
    /// \param[in] key Unique key for parameter
    /// \param[in] value Value to set the parameter
    template <typename T, typename U>
    inline ParametersList& setAs(const std::string& key, const U& value) {
      return set<T>(key, static_cast<T>(value));
    }
    ParametersList& rename(const std::string&, const std::string&);  ///< Rename the key to a parameter value
    ParametersList& operator+=(const ParametersList& oth);           ///< Concatenate two parameters containers
    ParametersList operator+(const ParametersList& oth) const;       ///< Concatenation of two parameters containers
    bool empty() const;                                              ///< Is the list empty?

    template <typename T>
    std::vector<std::string> keysOf() const;  ///< List of keys for one type in this list of parameters
    /// List of keys handled in this list of parameters
    /// \param[in] name_key Include the name variable?
    std::vector<std::string> keys(bool name_key = true) const;
    /// Get a string-converted version of a value
    /// \param[in] key Unique key for parameter
    /// \param[in] wrap Encapsulate the value with type()
    std::string getString(const std::string& key, bool wrap = false) const;
    /// Get a string-converted version of the module name if any
    /// \param[in] wrap Encapsulate the value with type()
    inline std::string getNameString(bool wrap = false) const { return getString(MODULE_NAME, wrap); }
    std::string serialise() const;  ///< Serialise a parameters collection into a parseable string

    friend std::ostream& operator<<(std::ostream&,
                                    const ParametersList&);  ///< Human-readable version of a parameters container
    const ParametersList& print(std::ostream&) const;        ///< Debugging-like printout of a parameters container
    std::string print(bool compact = false) const;           ///< Normal printout of a parameters container

  private:
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
  size_t ParametersList::erase<type>(const std::string&);
  REGISTER_CONTENT_TYPE
#undef __TYPE_ENUM
}  // namespace cepgen

#endif
