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

#ifndef CepGen_Utils_Message_h
#define CepGen_Utils_Message_h

#include <set>
#include <sstream>
#include <unordered_map>

#include "CepGen/Utils/Logger.h"

namespace cepgen {
  /// A generic message type
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date 27 Mar 2015
  struct Message {
    explicit Message() = default;  ///< Generic message constructor
    virtual ~Message() = default;
    /// Dump the full exception information in a given output stream
    /// \param[inout] os the output stream where the information is dumped
    virtual void dump(std::ostream* os = nullptr) const = 0;
    static std::string now();  ///< Human-readable date/time
  };

  /// A simple exception handler
  /// \date 24 Mar 2015
  class LoggedMessage : public Message {
  public:
    /// Enumeration of message type
    enum class MessageType {
      undefined = -1,  ///< Irregular message
      debug,           ///< Debugging information to be enabled
      verbatim,        ///< Raw information
      info,            ///< Prettified information
      warning,         ///< Casual non-stopping warning
    };
    /// Generic constructor
    /// \param[in] module exception classifier
    /// \param[in] from method invoking the exception
    /// \param[in] type exception type
    /// \param[in] file file where this occurred
    /// \param[in] lineno Line number where exception occurred
    explicit LoggedMessage(const std::string& module,
                           const std::string& from = "",
                           MessageType type = MessageType::undefined,
                           const std::string& file = "",
                           short lineno = 0) noexcept(true);
    LoggedMessage(const LoggedMessage&) noexcept(true);  ///< Copy constructor
    ~LoggedMessage() noexcept override;                  ///< Default destructor

    friend std::ostream& operator<<(std::ostream&, const MessageType&);  ///< Printout operator for message type

    //----- Overloaded stream operators

    /// Generic templated message feeder operator
    template <typename T>
    friend const LoggedMessage& operator<<(const LoggedMessage& log, const T& message) noexcept {
      auto& nc_except = const_cast<LoggedMessage&>(log);
      nc_except.message_ << message;
      return log;
    }
    /// Specialised feeder operator for booleans
    friend const LoggedMessage& operator<<(const LoggedMessage&, const bool&) noexcept;
    /// Specialised feeder operator for wide strings
    friend const LoggedMessage& operator<<(const LoggedMessage&, const std::wstring&) noexcept;
    /// Generic templated pair-variables feeder operator
    template <typename T, typename U>
    friend const LoggedMessage& operator<<(const LoggedMessage& log, const std::pair<T, U>& pair_var) noexcept {
      return log << "(" << pair_var.first << ", " << pair_var.second << ")";
    }
    /// Generic templated vector-variables feeder operator
    template <typename T>
    friend const LoggedMessage& operator<<(const LoggedMessage& log, const std::set<T>& set_var) noexcept {
      (void)(log << "[");
      std::string sep;
      if (!set_var.empty())
        for (const auto& var : set_var)
          log << sep << var, sep = ", ";
      return log << "]";
    }
    /// Generic templated vector-variables feeder operator
    template <typename T>
    friend const LoggedMessage& operator<<(const LoggedMessage& log, const std::vector<T>& vec_var) noexcept {
      (void)(log << "{");
      std::string sep;
      if (!vec_var.empty())
        for (const auto& var : vec_var)
          log << sep << var, sep = ", ";
      return log << "}";
    }
    /// Generic templated vector-variables feeder operator
    template <typename T, std::size_t N>
    friend const LoggedMessage& operator<<(const LoggedMessage& log, const std::array<T, N>& vec_var) noexcept {
      (void)(log << "{");
      std::string sep;
      if (!vec_var.empty())
        for (const auto& var : vec_var)
          log << sep << var, sep = ", ";
      return log << "}";
    }
    /// Generic templated mapping-variables feeder operator
    template <typename T, typename U>
    friend const LoggedMessage& operator<<(const LoggedMessage& log, const std::map<T, U>& map_var) noexcept {
      (void)(log << "{");
      std::string sep;
      if (!map_var.empty())
        for (const auto& var : map_var)
          log << sep << "{" << var.first << " -> " << var.second << "}", sep = ", ";
      return log << "}";
    }
    /// Generic templated mapping-variables feeder operator
    template <typename T, typename U>
    friend const LoggedMessage& operator<<(const LoggedMessage& log, const std::unordered_map<T, U>& map_var) noexcept {
      (void)(log << "{");
      std::string sep;
      if (!map_var.empty())
        for (const auto& var : map_var)
          log << sep << "{" << var.first << " -> " << var.second << "}", sep = ", ";
      return log << "}";
    }
    /// Pipe modifier operator
    friend const LoggedMessage& operator<<(const LoggedMessage& log, std::ios_base& (*f)(std::ios_base&)) noexcept {
      auto& nc_except = const_cast<LoggedMessage&>(log);
      f(nc_except.message_);
      return log;
    }

    /// Lambda function handler
    template <typename T>
    const LoggedMessage& log(T&& lambda) noexcept {
      lambda(*this);
      return *this;
    }

    std::string message() const { return message_.str(); }          ///< Human-readable message
    const std::string& from() const { return from_; }               ///< Origin of the exception
    const std::string& file() const { return file_; }               ///< File where the exception occurred
    short lineNumber() const { return line_num_; }                  ///< Line number where the exception occurred
    const MessageType& type() const { return type_; }               ///< Message type
    void dump(std::ostream* os = nullptr) const noexcept override;  ///< Human-readable dump of the message
    std::ostream& stream() { return message_; }                     ///< Output stream object

  protected:
    std::ostringstream message_;  ///< Message to throw
    std::string from_;            ///< Origin of the exception
    std::string file_;            ///< File
    short line_num_;              ///< Line number

  private:
    MessageType type_;    ///< Message type
    std::string module_;  ///< Exception classification
  };

  /// Placeholder for debugging messages if logging threshold is not reached
  /// \date Apr 2018
  struct NullStream : Message {
    using Message::Message;
    NullStream() = default;              ///< Empty constructor
    NullStream(const LoggedMessage&) {}  ///< Empty constructor
    void dump(std::ostream* = nullptr) const override {}
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
  };
}  // namespace cepgen

#ifdef _WIN32
#define __FUNC__ __FUNCSIG__
#else
#define __FUNC__ __PRETTY_FUNCTION__
#endif

#define CG_LOG                                                                    \
  (cepgen::utils::Logger::get().level() <= cepgen::utils::Logger::Level::nothing) \
      ? cepgen::NullStream()                                                      \
      : cepgen::LoggedMessage("Logging", __FUNC__, cepgen::LoggedMessage::MessageType::verbatim, __FILE__, __LINE__)
#define CG_INFO(mod)                \
  (!CG_LOG_MATCH(mod, information)) \
      ? cepgen::NullStream()        \
      : cepgen::LoggedMessage(mod, __FUNC__, cepgen::LoggedMessage::MessageType::info, __FILE__, __LINE__)
#define CG_DEBUG(mod)         \
  (!CG_LOG_MATCH(mod, debug)) \
      ? cepgen::NullStream()  \
      : cepgen::LoggedMessage(mod, __FUNC__, cepgen::LoggedMessage::MessageType::debug, __FILE__, __LINE__)
#define CG_DEBUG_LOOP(mod)              \
  (!CG_LOG_MATCH(mod, debugInsideLoop)) \
      ? cepgen::NullStream()            \
      : cepgen::LoggedMessage(mod, __FUNC__, cepgen::LoggedMessage::MessageType::debug, __FILE__, __LINE__)
#define CG_WARNING(mod)         \
  (!CG_LOG_MATCH(mod, warning)) \
      ? cepgen::NullStream()    \
      : cepgen::LoggedMessage(mod, __FUNC__, cepgen::LoggedMessage::MessageType::warning, __FILE__, __LINE__)

#endif
