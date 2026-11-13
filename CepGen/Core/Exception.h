/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2026  Laurent Forthomme
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

#include <memory>
#include <stdexcept>
#include <type_traits>

#include "CepGen/Utils/Message.h"

namespace cepgen {
  /// Standard exception message
  class Exception : public Message, public std::runtime_error {
  public:
    /// Enumeration of exception severities
    enum struct Type {
      undefined = -1,  ///< Irregular message
      error,           ///< General non-stopping error
      fatal            ///< Critical and stopping error
    };
    explicit Exception(const std::string& module_name,
                       const std::string& from = "",
                       Type type = Type::undefined,
                       const std::string& file = "",
                       int lineno = 0);
    Exception(const Exception&) noexcept = default;             ///< Copy constructor (sharing the state)
    Exception& operator=(const Exception&) noexcept = default;  ///< Assignment operator (sharing the state)
    ~Exception() noexcept override;                             ///< Destructor

    friend std::ostream& operator<<(std::ostream&, const Type&);  ///< Printout operator for the exception type

    /// Generic message feeder operator
    template <typename T>
    friend const Exception& operator<<(const Exception& exception, const T& message) {
      exception.messageBuffer() << message;
      return exception;
    }
    /// Lambda function handler, receiving a LoggedMessage reference to be fed
    template <typename F>
    const Exception& log(F&& lambda) const {
      lambda(messageBuffer());
      return *this;
    }

    /// Human-readable dump of the exception
    /// \param[inout] os Output stream where to dump the information. If not specified, the logger output is used,
    ///   and the exception is only written to it once, whatever the number of dump() calls
    void dump(std::ostream* os = nullptr) const noexcept override;
    /// Mark the exception as handled, preventing any automatic printout
    /// \note A fatal exception is always reported once, as it stops the process
    void handled() const noexcept;
    bool reported() const noexcept;  ///< Was this exception already reported (dumped or marked as handled)?
    const char* what() const noexcept override;

    std::string message() const;               ///< Human-readable message
    const std::string& from() const noexcept;  ///< Origin of the exception
    const std::string& file() const noexcept;  ///< File where the exception occurred
    int lineNumber() const noexcept;           ///< Line number where the exception occurred
    Type type() const noexcept;                ///< Exception severity
    std::ostream& stream() const;              ///< Output stream object

    /// Action triggered when the handling of a fatal exception is over
    using FatalHandler = void (*)();
    /// Replace the action triggered once a fatal exception is over (std::exit(EXIT_FAILURE) by default)
    /// \param[in] handler New action, or nullptr to restore the default one. If it returns, the execution continues
    ///   (e.g. for tests, or when CepGen is embedded in an application handling its own shutdown)
    /// \return Previous action
    static FatalHandler setFatalHandler(FatalHandler handler) noexcept;

  private:
    LoggedMessage& messageBuffer() const noexcept;  ///< Silent message buffer shared between all copies

    struct State;
    std::shared_ptr<State> state_;
  };
  static_assert(std::is_nothrow_copy_constructible_v<Exception>, "Exception must be nothrow copy-constructible");
}  // namespace cepgen

#define CG_ERROR(mod) cepgen::Exception(mod, __FUNC__, cepgen::Exception::Type::error, __FILE__, __LINE__)
#define CG_FATAL(mod) cepgen::Exception(mod, __FUNC__, cepgen::Exception::Type::fatal, __FILE__, __LINE__)
#define CG_ASSERT(assertion) \
  if (!(assertion))          \
    throw CG_FATAL("Assertion") << "Assertion '" << #assertion << "' failed.";

#endif
