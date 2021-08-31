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

#include <csignal>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  LoggedException::LoggedException(const char* module, Type type, const char* file, short lineno)
      : module_(module), type_(type), file_(file), line_num_(lineno) {}

  LoggedException::LoggedException(const char* from, const char* module, Type type, const char* file, short lineno)
      : from_(from), module_(module), type_(type), file_(file), line_num_(lineno) {}

  LoggedException::LoggedException(const LoggedException& rhs) noexcept
      : from_(rhs.from_),
        module_(rhs.module_),
        message_(rhs.message_.str()),
        type_(rhs.type_),
        line_num_(rhs.line_num_) {}

  LoggedException::~LoggedException() noexcept {
    if (type_ != Type::undefined)
      dump();
    // we stop this process' execution on fatal exception
    if (type_ == Type::fatal && raise(SIGINT) != 0)
      exit(0);
  }

  const LoggedException& operator<<(const LoggedException& exc, const std::wstring& var) {
    LoggedException& nc_except = const_cast<LoggedException&>(exc);
    nc_except.message_ << utils::tostring(var);
    return exc;
  }

  const char* LoggedException::what() const noexcept {
    (*utils::Logger::get().output) << "\n" << message_.str() << "\n";
    return from_.c_str();
  }

  std::ostream& LoggedException::dump(std::ostream& os) const {
    if (!utils::Logger::get().output)
      return os;

    switch (type_) {
      case Type::info:
        return os << type_ << ":\t" << message_.str() << "\n";
      case Type::debug:
        return os << type_ << " "
                  << utils::colourise(
                         from_, utils::Colour::yellow, utils::Modifier::underline | utils::Modifier::dimmed)
                  << (utils::Logger::get().extended()
                          ? " " +
                                utils::colourise(
                                    file_,
                                    utils::Colour::reset,
                                    utils::Modifier::bold | utils::Modifier::italic | utils::Modifier::dimmed) +
                                " @" +
                                utils::colourise(std::to_string(line_num_),
                                                 utils::Colour::reset,
                                                 utils::Modifier::italic | utils::Modifier::dimmed) +
                                "\n"
                          : ": ")
                  << utils::colourise(message_.str(), utils::Colour::reset, utils::Modifier::dimmed) << "\n";
      case Type::warning:
        return os << type_ << " "
                  << utils::colourise(from_, utils::Colour::reset, utils::Modifier::underline | utils::Modifier::dimmed)
                  << (utils::Logger::get().extended()
                          ? " " +
                                utils::colourise(
                                    file_,
                                    utils::Colour::reset,
                                    utils::Modifier::bold | utils::Modifier::italic | utils::Modifier::dimmed) +
                                " @" +
                                utils::colourise(std::to_string(line_num_),
                                                 utils::Colour::reset,
                                                 utils::Modifier::italic | utils::Modifier::dimmed)
                          : "")
                  << "\n\t" << message_.str() << "\n";
      case Type::verbatim:
        return os << message_.str() << "\n";
      case Type::undefined:
      case Type::error:
      case Type::fatal: {
        const std::string sep(80, '-');
        os << sep << "\n" << type_ << " occured at " << now() << "\n";
        if (!from_.empty())
          os << "  raised by: " << utils::colourise(from_, utils::Colour::reset, utils::Modifier::underline) << "\n";
        if (utils::Logger::get().extended()) {
          os << "  file: " << utils::colourise(file_, utils::Colour::reset, utils::Modifier::dimmed) << "\n";
          if (line_num_ != 0)
            os << "  line #" << line_num_ << "\n";
        }
        os << "\n" << message_.str() << "\n";
        return os << sep << "\n";
      }
    }
    return os;
  }

  char* LoggedException::now() {
    static char buffer[10];
    time_t rawtime;
    time(&rawtime);
    struct tm* timeinfo = localtime(&rawtime);
    strftime(buffer, 10, "%H:%M:%S", timeinfo);
    return buffer;
  }

  std::ostream& operator<<(std::ostream& os, const Exception::Type& type) {
    switch (type) {
      case Exception::Type::info:
        return os << utils::colourise("Info", utils::Colour::green, utils::Modifier::bold);
      case Exception::Type::debug:
        return os << utils::colourise("Debug", utils::Colour::yellow, utils::Modifier::bold);
      case Exception::Type::warning:
        return os << utils::colourise("Warning", utils::Colour::blue, utils::Modifier::bold);
      case Exception::Type::verbatim:
        return os << utils::colourise("Verbatim", utils::Colour::reset, utils::Modifier::bold);
      case Exception::Type::undefined:
        return os << utils::colourise("Undef'd exception", utils::Colour::reset, utils::Modifier::reverse);
      case Exception::Type::error:
        return os << utils::colourise("Error", utils::Colour::red, utils::Modifier::bold);
      case Exception::Type::fatal:
        return os << utils::colourise("Fatal error", utils::Colour::red, utils::Modifier::bold);
    }
    return os;
  }
}  // namespace cepgen
