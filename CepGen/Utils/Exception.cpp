#include "CepGen/Core/Exception.h"

#include <csignal>

#include "CepGen/Utils/String.h"

namespace cepgen {
  LoggedException::LoggedException(const char* module, Type type, short id)
      : module_(module), type_(type), error_num_(id) {}

  LoggedException::LoggedException(const char* from, const char* module, Type type, short id)
      : from_(from), module_(module), type_(type), error_num_(id) {}

  LoggedException::LoggedException(const LoggedException& rhs) noexcept
      : from_(rhs.from_),
        module_(rhs.module_),
        message_(rhs.message_.str()),
        type_(rhs.type_),
        error_num_(rhs.error_num_) {}

  LoggedException::~LoggedException() noexcept {
    if (type_ != Type::undefined)
      dump();
    // we stop this process' execution on fatal exception
    if (type_ == Type::fatal && raise(SIGINT) != 0)
      exit(0);
  }

  const char* LoggedException::what() const noexcept {
    (*utils::Logger::get().output) << "\n" << message_.str() << "\n";
    return from_.c_str();
  }

  std::string LoggedException::message() const { return message_.str(); }

  std::string LoggedException::from() const { return from_; }

  int LoggedException::errorNumber() const { return error_num_; }

  Exception::Type LoggedException::type() const { return type_; }

  std::ostream& LoggedException::dump(std::ostream& os) const {
    if (!utils::Logger::get().output)
      return os;

    switch (type_) {
      case Type::info:
        return os << type_ << ":\t" << message_.str() << "\n";
      case Type::debug:
      case Type::warning:
        return os << type_ << " " << utils::colourise(from_, utils::Colour::reset, utils::Modifier::underline) << "\n\t"
                  << message_.str() << "\n";
      case Type::verbatim:
        return os << message_.str() << "\n";
      case Type::undefined:
      case Type::error:
      case Type::fatal: {
        const std::string sep(80, '-');
        os << sep << "\n" << type_ << " occured at " << now() << "\n";
        if (!from_.empty())
          os << "  raised by: " << utils::colourise(from_, utils::Colour::reset, utils::Modifier::underline) << "\n";
        if (errorNumber() != 0)
          os << "  error #" << error_num_ << "\n";
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
        return os << utils::colourise("Debug", utils::Colour::yellow, utils::Modifier::reverse);
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
