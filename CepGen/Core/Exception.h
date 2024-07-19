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
  class Exception : public LoggedMessage, public std::runtime_error {
  public:
    /// Enumeration of exception severities
    enum Type {
      undefined = -1,  ///< Irregular message
      error,           ///< General non-stopping error
      fatal            ///< Critical and stopping error
    };
    explicit Exception(const char* mod,
                       const char* from = "",
                       Type type = Type::undefined,
                       const char* file = "",
                       short lineno = 0) noexcept;
    /// Destructor (potentially killing the process)
    virtual ~Exception() noexcept override;

    /// Printout operator for exception type
    friend std::ostream& operator<<(std::ostream&, const Type&);

    template <typename T>
    inline friend const Exception& operator<<(const Exception& exc, const T& var) noexcept {
      (const LoggedMessage&)exc << var;
      return exc;
    }
    /// Human-readable dump of the exception
    void dump(std::ostream* os = nullptr) const noexcept override;

    const char* what() const noexcept override;

  private:
    Type type_;
    mutable char what_[50];
  };
}  // namespace cepgen

#define CG_ERROR(mod) cepgen::Exception(mod, __FUNC__, cepgen::Exception::Type::error, __FILE__, __LINE__)
#define CG_FATAL(mod) cepgen::Exception(mod, __FUNC__, cepgen::Exception::Type::fatal, __FILE__, __LINE__)
#define CG_ASSERT(assertion) \
  if (!(assertion))          \
    throw CG_FATAL("Assertion") << "Assertion '" << #assertion << " failed.";

#endif
