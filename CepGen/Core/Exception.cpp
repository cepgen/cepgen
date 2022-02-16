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

#include <csignal>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  Exception::Exception(const char* mod, const char* from, Type type, const char* file, short lineno) noexcept
      : LoggedMessage(mod, from, MessageType::undefined, file, lineno),
        std::runtime_error("cepgen::Exception"),
        type_(type) {}

  Exception::Exception(const Exception& oth) noexcept
      : LoggedMessage(oth), std::runtime_error(oth.what()), type_(oth.type_) {}

  Exception::~Exception() noexcept {
    dump();
    // we stop this process' execution on fatal exception
    if (type_ == Type::fatal && raise(SIGINT) != 0)
      exit(0);
  }

  const char* Exception::what() const noexcept { return message_.str().c_str(); }

  void Exception::dump(std::ostream* os) const noexcept {
    if (!os)
      os = utils::Logger::get().output;
    if (!os)
      return;

    const std::string sep(80, '-');
    (*os) << sep << "\n" << type_ << " occured at " << now() << "\n";
    if (!from_.empty())
      (*os) << "  raised by: " << utils::colourise(from_, utils::Colour::none, utils::Modifier::underline) << "\n";
    if (utils::Logger::get().extended() && !file_.empty()) {
      (*os) << "  file: " << utils::colourise(file_, utils::Colour::none, utils::Modifier::dimmed) << "\n";
      if (line_num_ != 0)
        (*os) << "  line #" << line_num_ << "\n";
    }
    (*os) << "\n" << message_.str() << "\n" << sep << "\n";
  }

  std::ostream& operator<<(std::ostream& os, const Exception::Type& type) {
    switch (type) {
      case Exception::Type::error:
        return os << utils::colourise("Error", utils::Colour::red, utils::Modifier::bold);
      case Exception::Type::fatal:
        return os << utils::colourise("Fatal error", utils::Colour::red, utils::Modifier::bold);
      case Exception::Type::undefined:
        return os << utils::colourise("Undef'd exception", utils::Colour::none, utils::Modifier::reverse);
    }
    return os;
  }
}  // namespace cepgen
