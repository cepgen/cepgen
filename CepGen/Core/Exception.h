/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#ifndef CepGen_Core_Exception_h
#define CepGen_Core_Exception_h

#include "CepGen/Utils/Message.h"

namespace cepgen {
  class Exception final : public LoggedMessage, public std::runtime_error {
  public:
    /// Enumeration of exception severities
    enum Type {
      undefined = -1,  ///< Irregular message
      error,           ///< General non-stopping error
      fatal            ///< Critical and stopping error
    };
    explicit Exception(
        const char* mod, const char* from = "", Type type = Type::undefined, const char* file = "", short lineno = 0);
    Exception(Exception&&);
    /// Destructor (potentially killing the process)
    virtual ~Exception() override;

    /// Printout operator for exception type
    friend std::ostream& operator<<(std::ostream&, const Type&);
    /// Dump the full exception information in a given output stream
    /// \param[inout] os the output stream where the information is dumped
    virtual std::ostream& dump(std::ostream& os = *utils::Logger::get().output) const = 0;
    /// Exception message
    virtual std::string message() const = 0;
  };

  /// A simple exception handler
  /// \date 24 Mar 2015
  class LoggedException : public Exception, public std::exception {
  public:
    /// Generic constructor
    /// \param[in] module exception classifier
    /// \param[in] type exception type
    /// \param[in] lineno Line number where exception occured
    explicit LoggedException(const char* module = "",
                             Type type = Type::undefined,
                             const char* file = "",
                             short lineno = 0);
    /// Generic constructor
    /// \param[in] from method invoking the exception
    /// \param[in] module exception classifier
    /// \param[in] type exception type
    /// \param[in] lineno Line number where exception occured
    explicit LoggedException(
        const char* from, const char* module, Type type = Type::undefined, const char* file = "", short lineno = 0);
    /// Copy constructor
    LoggedException(const LoggedException& rhs) noexcept;
    /// Default destructor (potentially killing the process)
    ~LoggedException() noexcept override;

    //----- Overloaded stream operators

    /// Generic templated message feeder operator
    template <typename T>
    inline friend const LoggedException& operator<<(const LoggedException& exc, const T& var) {
      LoggedException& nc_except = const_cast<LoggedException&>(exc);
      nc_except.message_ << var;
      return exc;
    }
    /// Specialised feeder operator for booleans
    friend const LoggedException& operator<<(const LoggedException&, const bool&);
    /// Specialised feeder operator for wide strings
    friend const LoggedException& operator<<(const LoggedException&, const std::wstring&);
    /// Generic templated pair-variables feeder operator
    template <typename T, typename U>
    inline friend const LoggedException& operator<<(const LoggedException& exc, const std::pair<T, U>& pair_var) {
      return exc << "(" << pair_var.first << ", " << pair_var.second << ")";
    }
    /// Generic templated vector-variables feeder operator
    template <typename T>
    inline friend const LoggedException& operator<<(const LoggedException& exc, const std::vector<T>& vec_var) {
      exc << "{";
      std::string sep;
      if (!vec_var.empty())
        for (const auto& var : vec_var)
          exc << sep << var, sep = ", ";
      return exc << "}";
    }
    /// Generic templated vector-variables feeder operator
    template <typename T, std::size_t N>
    inline friend const LoggedException& operator<<(const LoggedException& exc, const std::array<T, N>& vec_var) {
      exc << "{";
      std::string sep;
      if (!vec_var.empty())
        for (const auto& var : vec_var)
          exc << sep << var, sep = ", ";
      return exc << "}";
    }
    /// Generic templated mapping-variables feeder operator
    template <typename T, typename U>
    inline friend const LoggedException& operator<<(const LoggedException& exc, const std::map<T, U>& map_var) {
      exc << "{";
      std::string sep;
      if (!map_var.empty())
        for (const auto& var : map_var)
          exc << sep << "{" << var.first << " -> " << var.second << "}", sep = ", ";
      return exc << "}";
    }
    /// Pipe modifier operator
    inline friend const LoggedException& operator<<(const LoggedException& exc, std::ios_base& (*f)(std::ios_base&)) {
      LoggedException& nc_except = const_cast<LoggedException&>(exc);
      f(nc_except.message_);
      return exc;
    }

    /// Lambda function handler
    template <typename T>
    inline LoggedException& log(T&& lam) {
      lam(*this);
      return *this;
    }

    const char* what() const noexcept override;

  private:
    Type type_;
  };
}  // namespace cepgen

#define CG_ERROR(mod) cepgen::Exception(mod, __FUNC__, cepgen::Exception::Type::error, __FILE__, __LINE__)
#define CG_FATAL(mod) cepgen::Exception(mod, __FUNC__, cepgen::Exception::Type::fatal, __FILE__, __LINE__)

#endif
