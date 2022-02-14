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

#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  char* Message::now() {
    static char buffer[10];
    time_t rawtime;
    time(&rawtime);
    struct tm* timeinfo = localtime(&rawtime);
    strftime(buffer, 10, "%H:%M:%S", timeinfo);
    return buffer;
  }

  LoggedMessage::LoggedMessage(
      const char* mod, const char* from, MessageType type, const char* file, short lineno) noexcept
      : from_(from), file_(file), line_num_(lineno), type_(type), module_(mod) {}

  LoggedMessage::LoggedMessage(const LoggedMessage& rhs) noexcept
      : message_(rhs.message_.str()),  // only reason to customise the copy constructor
        from_(rhs.from_),
        file_(rhs.file_),
        line_num_(rhs.line_num_),
        type_(rhs.type_),
        module_(rhs.module_) {}

  LoggedMessage::~LoggedMessage() noexcept {
    if (type_ != MessageType::undefined)
      dump();
  }

  const LoggedMessage& operator<<(const LoggedMessage& exc, const bool& var) {
    LoggedMessage& nc_except = const_cast<LoggedMessage&>(exc);
    nc_except.message_ << (var ? utils::colourise("true", utils::Colour::green)
                               : utils::colourise("false", utils::Colour::red));
    return exc;
  }

  const LoggedMessage& operator<<(const LoggedMessage& exc, const std::wstring& var) {
    LoggedMessage& nc_except = const_cast<LoggedMessage&>(exc);
    nc_except.message_ << utils::tostring(var);
    return exc;
  }

  void LoggedMessage::dump(std::ostream* os) const {
    if (!os)
      os = utils::Logger::get().output;
    if (!os)
      return;

    if (type_ == MessageType::verbatim || type_ == MessageType::undefined) {
      (*os) << message_.str() << "\n";
      return;
    }
    (*os) << type_;
    if (type_ == MessageType::info) {
      if (utils::Logger::get().extended())
        (*os) << utils::colourise(
            "{" + from_ + "}\n\t", utils::Colour::none, utils::Modifier::dimmed | utils::Modifier::italic);
      else
        (*os) << ":\t";
      (*os) << message_.str() << "\n";
      return;
    }
    (*os) << " " << utils::colourise(from_, utils::Colour::none, utils::Modifier::underline | utils::Modifier::dimmed);
    if (utils::Logger::get().extended())
      (*os) << " "
            << utils::colourise(file_,
                                utils::Colour::none,
                                utils::Modifier::bold | utils::Modifier::italic | utils::Modifier::dimmed)
            << " @"
            << utils::colourise(
                   std::to_string(line_num_), utils::Colour::none, utils::Modifier::italic | utils::Modifier::dimmed)
            << "\n";
    else
      (*os) << ": ";
    if (type_ == MessageType::debug) {
      (*os) << utils::colourise(message_.str(), utils::Colour::none, utils::Modifier::dimmed) << "\n";
      return;
    }
    if (type_ == MessageType::warning) {
      (*os) << message_.str() << "\n";
      return;
    }
  }

  std::ostream& operator<<(std::ostream& os, const LoggedMessage::MessageType& type) {
    switch (type) {
      case LoggedMessage::MessageType::info:
        return os << utils::colourise("Info", utils::Colour::green, utils::Modifier::bold);
      case LoggedMessage::MessageType::debug:
        return os << utils::colourise("Debug", utils::Colour::yellow, utils::Modifier::bold);
      case LoggedMessage::MessageType::warning:
        return os << utils::colourise("Warning", utils::Colour::blue, utils::Modifier::bold);
      case LoggedMessage::MessageType::verbatim:
        return os << utils::colourise("Verbatim", utils::Colour::none, utils::Modifier::bold);
      case LoggedMessage::MessageType::undefined:
        return os << utils::colourise("Undef'd exception", utils::Colour::none, utils::Modifier::reverse);
    }
    return os;
  }
}  // namespace cepgen
