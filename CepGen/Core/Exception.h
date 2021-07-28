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
    /// \param[in] id exception code (useful for logging)
    explicit LoggedException(const char* module = "", Type type = Type::undefined, short id = 0);
    /// Generic constructor
    /// \param[in] from method invoking the exception
    /// \param[in] module exception classifier
    /// \param[in] type exception type
    /// \param[in] id exception code (useful for logging)
    explicit LoggedException(const char* from, const char* module, Type type = Type::undefined, short id = 0);
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
    /// Specialised feeder operator for wide strings
    friend const LoggedException& operator<<(const LoggedException& exc, const std::wstring& var);
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
    std::string message() const override;

    /// Origin of the exception
    std::string from() const;
    /// Exception code
    int errorNumber() const;
    /// Exception type
    Type type() const;
    /// Human-readable dump of the exception message
    std::ostream& dump(std::ostream& os = *utils::Logger::get().output) const override;

  private:
    static char* now();
    /// Origin of the exception
    std::string from_;
    /// Exception classification
    std::string module_;
    /// Message to throw
    std::ostringstream message_;
    /// Exception type
    Type type_;
    /// Integer exception number
    short error_num_;
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

#define CG_LOG cepgen::LoggedException(__FUNC__, "Logging", cepgen::Exception::Type::verbatim)
#define CG_INFO(mod)                                       \
  (!CG_LOG_MATCH(mod, information)) ? cepgen::NullStream() \
                                    : cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::info)
#define CG_DEBUG(mod)                                \
  (!CG_LOG_MATCH(mod, debug)) ? cepgen::NullStream() \
                              : cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::debug)
#define CG_DEBUG_LOOP(mod)                                     \
  (!CG_LOG_MATCH(mod, debugInsideLoop)) ? cepgen::NullStream() \
                                        : cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::debug)
#define CG_WARNING(mod)                                \
  (!CG_LOG_MATCH(mod, warning)) ? cepgen::NullStream() \
                                : cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::warning)
#define CG_ERROR(mod) cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::error)
#define CG_FATAL(mod) cepgen::LoggedException(__FUNC__, mod, cepgen::Exception::Type::fatal)

#endif
