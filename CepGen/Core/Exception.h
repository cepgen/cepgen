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

#ifndef CepGen_Core_Exception_h
#define CepGen_Core_Exception_h

#include "CepGen/Utils/Logger.h"

namespace cepgen {
  /// A generic exception type
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date 27 Mar 2015
  struct Exception {
    /// Generic exception constructor
    explicit inline Exception() = default;
    virtual ~Exception() noexcept = default;
    /// Enumeration of exception severities
    enum class Type {
      undefined = -1,  ///< Irregular exception
      debug,           ///< Debugging information to be enabled
      verbatim,        ///< Raw information
      info,            ///< Prettified information
      warning,         ///< Casual non-stopping warning
      error,           ///< General non-stopping error
      fatal            ///< Critical and stopping error
    };
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
  class LoggedException final : public Exception, public std::exception {
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
    std::string message() const override { return message_.str(); }

    /// Origin of the exception
    const std::string& from() const { return from_; }
    /// File where the exception occured
    const std::string& file() const { return file_; }
    /// Line number where the exception occured
    short lineNumber() const { return line_num_; }
    /// Exception type
    const Type& type() const { return type_; }
    /// Human-readable dump of the exception message
    std::ostream& dump(std::ostream& os = *utils::Logger::get().output) const override;

  private:
    static char* now();
    std::string from_;            ///< Origin of the exception
    std::string module_;          ///< Exception classification
    std::ostringstream message_;  ///< Message to throw
    Type type_;                   ///< Exception type
    std::string file_;            ///< File
    short line_num_;              ///< Line number
  };

  /// Placeholder for debugging messages if logging threshold is not reached
  /// \date Apr 2018
  struct NullStream : Exception {
    using Exception::Exception;
    /// Empty constructor
    inline NullStream() {}
    /// Empty constructor
    inline NullStream(const LoggedException&) {}
    std::ostream& dump(std::ostream& os) const override { return os; }
    /// Stream operator (null and void)
    template <class T>
    NullStream& operator<<(const T&) {
      return *this;
    }
    /// Lambda function handler (null and void)
    template <typename T>
    NullStream& log(T&&) {
      return *this;
    }
    std::string message() const override { return ""; }
  };
}  // namespace cepgen

#ifdef _WIN32
#define __FUNC__ __FUNCSIG__
#else
#define __FUNC__ __PRETTY_FUNCTION__
#endif

#define CG_LOG cepgen::LoggedException(__FUNC__, "Logging", cepgen::Exception::Type::verbatim, __FILE__, __LINE__)
#define CG_INFO(mod)                \
  (!CG_LOG_MATCH(mod, information)) \
      ? cepgen::NullStream()        \
      : cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::info, __FILE__, __LINE__)
#define CG_DEBUG(mod)         \
  (!CG_LOG_MATCH(mod, debug)) \
      ? cepgen::NullStream()  \
      : cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::debug, __FILE__, __LINE__)
#define CG_DEBUG_LOOP(mod)              \
  (!CG_LOG_MATCH(mod, debugInsideLoop)) \
      ? cepgen::NullStream()            \
      : cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::debug, __FILE__, __LINE__)
#define CG_WARNING(mod)         \
  (!CG_LOG_MATCH(mod, warning)) \
      ? cepgen::NullStream()    \
      : cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::warning, __FILE__, __LINE__)
#define CG_ERROR(mod) cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::error, __FILE__, __LINE__)
#define CG_FATAL(mod) cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::fatal, __FILE__, __LINE__)

#endif
