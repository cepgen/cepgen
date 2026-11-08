/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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
  /// Standard exception message
  class Exception : public LoggedMessage, public std::runtime_error {
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
                       short lineno = 0) noexcept(true);
    Exception(const Exception&) noexcept(true);  ///< Copy constructor
    ~Exception() noexcept override;              ///< Destructor (potentially killing the process)

    friend std::ostream& operator<<(std::ostream&, const Type&);  ///< Printout operator for the exception type

    template <typename T>
    friend const Exception& operator<<(const Exception& exception, const T& message) noexcept {
      static_cast<const LoggedMessage&>(exception) << message;
      return exception;
    }
    void dump(std::ostream* = nullptr) const noexcept(true) override;  ///< Human-readable dump of the exception
    const char* what() const noexcept(true) override;

  private:
    const Type type_{Type::undefined};
    mutable char what_[50]{""};
  };
  static_assert(std::is_nothrow_copy_constructible_v<Exception>, "Exception must be nothrow copy-constructible");
}  // namespace cepgen

#define CG_ERROR(mod) cepgen::Exception(mod, __FUNC__, cepgen::Exception::Type::error, __FILE__, __LINE__)
#define CG_FATAL(mod) cepgen::Exception(mod, __FUNC__, cepgen::Exception::Type::fatal, __FILE__, __LINE__)
#define CG_ASSERT(assertion) \
  if (!(assertion))          \
    throw CG_FATAL("Assertion") << "Assertion '" << #assertion << "' failed.";

#endif
