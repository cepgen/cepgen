/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2021  Laurent Forthomme
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
#include <vector>

#include "CepGen/Utils/Limits.h"

#define DEFINE_TYPE(type)                                                        \
  template <>                                                                    \
  bool ParametersList::has<type>(const std::string& key) const;                  \
  template <>                                                                    \
  type ParametersList::get<type>(const std::string& key, const type& def) const; \
  template <>                                                                    \
  type& ParametersList::operator[]<type>(const std::string& key);                \
  template <>                                                                    \
  ParametersList& ParametersList::set<type>(std::string key, const type&);       \
  template <>                                                                    \
  std::vector<std::string> ParametersList::keysOf<type>() const;

namespace cepgen {
  /// Parameters container
  class ParametersList {
  private:
    /// Retrieve the default argument for a given variable type
    template <typename T>
    struct default_arg {
      /// Default variable argument
      static T get() { return T(); }
    };

  public:
    ParametersList() = default;
    /// Copy constructor
    ParametersList(const ParametersList&);
    ~ParametersList() {}  // required for unique_ptr initialisation! avoids cleaning all individual objects
    ParametersList& operator=(const ParametersList&) = default;  ///< Assignment operator
    /// Equality operator
    bool operator==(const ParametersList&) const;
    /// Inequality operator
    bool operator!=(const ParametersList& oth) const { return !operator==(oth); }
    /// Feed a control string to the list of parameters
    ParametersList& feed(const std::string&);
    /// Check if a given parameter is handled in this list
    template <typename T>
    bool has(const std::string& key) const;
    /// Erase a parameter with key
    /// \return Number of key-indexed values erased
    size_t erase(const std::string&);
    /// Does the parameters list have a name key?
    template <typename T>
    bool hasName() const {
      return has<T>(MODULE_NAME);
    }
    /// Retrieve the module name if any
    template <typename T>
    inline T name(const T& def = default_arg<T>::get()) const {
      if (!has<T>(MODULE_NAME))
        return def;
      return get<T>(MODULE_NAME);
    }
    /// Set the module name
    template <typename T>
    inline ParametersList& setName(const T& value) {
      return set<T>(MODULE_NAME, value);
    }
    /// Fill a variable with the key content if exists
    template <typename T>
    const ParametersList& fill(const std::string& key, T& value) const {
      if (has<T>(key))
        value = get<T>(key);
      return *this;
    }
    /// Get a parameter value
    template <typename T>
    T get(const std::string& key, const T& def = default_arg<T>::get()) const;
    /// Get a recasted parameter value
    template <typename T, typename U>
    inline U getAs(const std::string& key, const U& def = default_arg<U>::get()) const {
      return static_cast<U>(get<T>(key, static_cast<T>(def)));
    }
    /// Reference to a parameter value
    template <typename T>
    T& operator[](const std::string& key);
    /// Set a parameter value
    template <typename T>
    ParametersList& set(std::string key, const T& value);
    /// Set a recasted parameter value
    template <typename T, typename U>
    inline ParametersList& setAs(std::string key, const U& value) {
      return set<T>(std::move(key), static_cast<T>(value));
    }
    /// Rename the key to a parameter value
    ParametersList& rename(const std::string& old_key, const std::string& new_key);
    /// Concatenate two parameters containers
    ParametersList& operator+=(const ParametersList& oth);
    /// Concatenation of two parameters containers
    ParametersList operator+(const ParametersList& oth) const;
    /// Is the list empty?
    bool empty() const;

    /// List of keys for one type in this list of parameters
    template <typename T>
    std::vector<std::string> keysOf() const;
    /// List of keys handled in this list of parameters
    /// \param[in] name_key Include the name variable?
    std::vector<std::string> keys(bool name_key = true) const;
    /// Get a string-converted version of a value
    /// \param[in] wrap Encapsulate the value with type()
    std::string getString(const std::string& key, bool wrap = false) const;
    /// Serialise a parameters collection into a parseable string
    std::string serialise() const;

    /// Human-readable version of a parameters container
    friend std::ostream& operator<<(std::ostream& os, const ParametersList&);
    /// Debugging-like printout of a parameters container
    const ParametersList& print(std::ostream&) const;
    /// Indexing key for the module name
    static constexpr char MODULE_NAME[] = "mod_name";

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
    std::unordered_map<std::string, std::vector<ParametersList> > vec_param_values_;
    std::unordered_map<std::string, std::vector<std::vector<double> > > vec_vec_dbl_values_;
  };

  DEFINE_TYPE(ParametersList)
  DEFINE_TYPE(bool)
  DEFINE_TYPE(int)
  DEFINE_TYPE(unsigned long long)
  DEFINE_TYPE(double)
  DEFINE_TYPE(std::string)
  DEFINE_TYPE(Limits)
  DEFINE_TYPE(std::vector<int>)
  DEFINE_TYPE(std::vector<double>)
  DEFINE_TYPE(std::vector<std::string>)
  DEFINE_TYPE(std::vector<ParametersList>)
  DEFINE_TYPE(std::vector<std::vector<double> >)
}  // namespace cepgen

#undef DEFINE_TYPE
#endif
