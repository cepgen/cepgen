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

#include <atomic>
#include <cstdlib>
#include <exception>
#include <mutex>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

namespace {
  /// Default action once a fatal exception is over: stop the process
  [[noreturn]] void exitOnFatal() { std::exit(EXIT_FAILURE); }
  std::atomic<Exception::FatalHandler> fatal_handler{exitOnFatal};
}  // namespace

/// State shared by all copies of an exception
struct Exception::State {
  State(const std::string& module_name, const std::string& from, Type type, const std::string& file, int lineno) noexcept
      : message(module_name, from, LoggedMessage::MessageType::undefined, file, lineno), type(type) {}
  /// Last copy of the exception is gone: report it if nobody did it before, and stop the process if it is fatal
  ~State() {
    if (type >= Type::error && !reported.exchange(true))
      print(nullptr);
    if (type == Type::fatal)
      if (const auto handler = fatal_handler.load())
        handler();
  }
  State(const State&) = delete;
  State& operator=(const State&) = delete;

  void print(std::ostream* os) const noexcept {
    try {
      if (!os)
        os = utils::Logger::get().output().get();
      if (!os)
        return;
      static const std::string separator(80, '-');
      std::ostringstream buf;  // format the full message first, to write it at once
      buf << separator << "\n" << type << " occurred at " << Message::now() << "\n";
      if (!message.from().empty())
        buf << "  raised by: " << colourise(message.from(), utils::Colour::none, utils::Modifier::underline) << "\n";
      if (utils::Logger::get().extended() && !message.file().empty()) {
        buf << "  file: " << colourise(message.file(), utils::Colour::none, utils::Modifier::dimmed) << "\n";
        if (message.lineNumber() != 0)
          buf << "  line #" << message.lineNumber() << "\n";
      }
      buf << "\n" << message.message() << "\n" << separator << "\n";
      Message::write(*os, buf.str());
      os->flush();   // exceptions are rare, and may be followed by a process abortion: do not keep them buffered
    } catch (...) {  // never let a printout failure escape (this may be called from a destructor)
    }
  }

  LoggedMessage message;              ///< Silent (undefined type) message buffer, also holding the origin information
  const Type type;                    ///< Exception severity
  std::atomic<bool> reported{false};  ///< Was the exception already reported to the logger output?
  std::once_flag what_built;          ///< Has the what() string already been built?
  std::string what;                   ///< Cached what() string, alive as long as one copy of the exception is
};

Exception::Exception(
    const std::string& module_name, const std::string& from, Type type, const std::string& file, int lineno)
    : std::runtime_error("cepgen::Exception"), state_(std::make_shared<State>(module_name, from, type, file, lineno)) {}

Exception::~Exception() noexcept = default;

LoggedMessage& Exception::messageBuffer() const noexcept { return state_->message; }

const char* Exception::what() const noexcept {
  try {
    std::call_once(state_->what_built, [this] {
      state_->what = "cepgen::Exception from " + state_->message.from() + ":\n\t" + state_->message.message();
    });
    return state_->what.c_str();
  } catch (...) {
    return std::runtime_error::what();
  }
}

void Exception::dump(std::ostream* os) const noexcept {
  if (os) {  // explicit output stream requested: always write
    state_->reported = true;
    state_->print(os);
    return;
  }
  if (!state_->reported.exchange(true))  // logger output: only write once
    state_->print(nullptr);
}

void Exception::handled() const noexcept {
  if (state_->type != Type::fatal)  // fatal exception always reported before process ends
    state_->reported = true;
}

Exception::FatalHandler Exception::setFatalHandler(FatalHandler handler) noexcept {
  return fatal_handler.exchange(handler ? handler : exitOnFatal);
}

bool Exception::reported() const noexcept { return state_->reported; }

std::string Exception::message() const { return state_->message.message(); }

const std::string& Exception::from() const noexcept { return state_->message.from(); }

const std::string& Exception::file() const noexcept { return state_->message.file(); }

int Exception::lineNumber() const noexcept { return state_->message.lineNumber(); }

Exception::Type Exception::type() const noexcept { return state_->type; }

std::ostream& Exception::stream() const { return state_->message.stream(); }

namespace {
  /// Report uncaught exception before the process is terminated
  std::terminate_handler previous_terminate_handler = nullptr;
  [[noreturn]] void reportUncaughtException() noexcept {
    if (auto eptr = std::current_exception())
      try {
        std::rethrow_exception(eptr);
      } catch (const Exception& exception) {
        exception.dump();
      } catch (...) {
      }
    if (previous_terminate_handler)
      previous_terminate_handler();
    std::abort();
  }
  [[maybe_unused]] const bool terminate_handler_installed = [] {
    previous_terminate_handler = std::set_terminate(reportUncaughtException);
    return true;
  }();
}  // namespace

namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const Exception::Type& type) {
    switch (type) {
      case Exception::Type::error:
        return os << colourise("Error", utils::Colour::red, utils::Modifier::bold);
      case Exception::Type::fatal:
        return os << colourise("Fatal error", utils::Colour::red, utils::Modifier::bold);
      case Exception::Type::undefined:
        return os << colourise("Undefined exception", utils::Colour::none, utils::Modifier::reverse);
    }
    return os;
  }
}  // namespace cepgen
