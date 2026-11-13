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

#include <mutex>

#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

std::string Message::now() { return utils::timeAs("%H:%M:%S"); }

void Message::write(std::ostream& os, const std::string& formatted_message) noexcept {
  static std::mutex output_mutex;  // avoid interleaved printouts from concurrent threads
  try {
    const std::lock_guard<std::mutex> lock(output_mutex);
    os << formatted_message;
  } catch (...) {  // never let a printout failure escape (this may be called from a destructor)
  }
}

LoggedMessage::LoggedMessage(
    const std::string& module, const std::string& from, MessageType type, const std::string& file, int lineno) noexcept
    : from_(from), file_(file), line_num_(lineno), type_(type), module_(module) {}

LoggedMessage::~LoggedMessage() noexcept {
  if (type_ != MessageType::undefined)
    dump();
}

void LoggedMessage::dump(std::ostream* os) const noexcept {
  try {
    if (!os)
      os = utils::Logger::get().output().get();
    if (!os)
      return;

    std::ostringstream buf;  // format the full message first, to write it at once
    if (type_ == MessageType::verbatim || type_ == MessageType::undefined)
      buf << message_.str() << "\n";
    else {
      buf << type_;
      if (type_ == MessageType::info) {
        if (utils::Logger::get().extended())
          buf << colourise(
              "{" + from_ + "}\n\t", utils::Colour::none, utils::Modifier::dimmed | utils::Modifier::italic);
        else
          buf << ":\t";
        buf << message_.str() << "\n";
      } else {
        buf << " " << colourise(from_, utils::Colour::none, utils::Modifier::underline | utils::Modifier::dimmed);
        if (utils::Logger::get().extended())
          buf << " "
              << colourise(file_,
                           utils::Colour::none,
                           utils::Modifier::bold | utils::Modifier::italic | utils::Modifier::dimmed)
              << ":"
              << colourise(
                     std::to_string(line_num_), utils::Colour::none, utils::Modifier::italic | utils::Modifier::dimmed)
              << "\n";
        else
          buf << ": ";
        if (type_ == MessageType::debug)
          buf << colourise(message_.str(), utils::Colour::none, utils::Modifier::dimmed) << "\n";
        else if (type_ == MessageType::warning)
          buf << message_.str() << "\n";
      }
    }
    write(*os, buf.str());
  } catch (...) {  // never let a printout failure escape (this is called from the destructor)
  }
}

namespace cepgen {
  const LoggedMessage& operator<<(const LoggedMessage& exc, const bool& var) {
    LoggedMessage& nc_except = const_cast<LoggedMessage&>(exc);
    nc_except.message_ << (var ? colourise("true", utils::Colour::green) : colourise("false", utils::Colour::red));
    return exc;
  }

  const LoggedMessage& operator<<(const LoggedMessage& exc, const std::wstring& var) {
    LoggedMessage& nc_except = const_cast<LoggedMessage&>(exc);
    nc_except.message_ << utils::toString(var);
    return exc;
  }

  std::ostream& operator<<(std::ostream& os, const LoggedMessage::MessageType& type) {
    switch (type) {
      case LoggedMessage::MessageType::info:
        return os << colourise("Info", utils::Colour::green, utils::Modifier::bold);
      case LoggedMessage::MessageType::debug:
        return os << colourise("Debug", utils::Colour::yellow, utils::Modifier::bold);
      case LoggedMessage::MessageType::warning:
        return os << colourise("Warning", utils::Colour::blue, utils::Modifier::bold);
      case LoggedMessage::MessageType::verbatim:
        return os << colourise("Verbatim", utils::Colour::none, utils::Modifier::bold);
      case LoggedMessage::MessageType::undefined:
        return os << colourise("Undefined exception", utils::Colour::none, utils::Modifier::reverse);
    }
    return os;
  }
}  // namespace cepgen
