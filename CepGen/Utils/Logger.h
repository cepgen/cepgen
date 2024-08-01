/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2015-2024  Laurent Forthomme
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

#ifndef CepGen_Utils_Logger_h
#define CepGen_Utils_Logger_h

#include <regex>
#include <vector>

namespace cepgen::utils {
  /// General purpose message logger
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date 15 Oct 2015
  class Logger {
  public:
    static Logger& get(std::ostream* = nullptr);  ///< Retrieve the running instance of the logger

    /// \brief Add a new rule to display exceptions/messages
    /// \param[in] rule Regex rule to handle
    void addExceptionRule(const std::string& rule);
    /// Collection of logging exceptions
    inline const std::vector<std::regex>& exceptionRules() const { return allowed_exc_; }

    typedef std::unique_ptr<std::ostream, std::function<void(std::ostream*)> > StreamHandler;

    /// Logging threshold for the output stream
    enum class Level { nothing = 0, error, warning, information, debug, debugInsideLoop };
    /// Is the module set to be displayed/logged?
    /// \param[in] tmpl Module name to probe
    /// \param[in] lev Upper verbosity level
    bool passExceptionRule(const std::string& tmpl, const Level& lev) const;
    inline Level level() const { return level_; }                  ///< Logging threshold
    inline void setLevel(Level level) { level_ = level; }          ///< Set the logging threshold
    inline bool extended() const { return extended_; }             ///< Also show extended information?
    inline void setExtended(bool ext = true) { extended_ = ext; }  ///< Set the extended information flag
    bool isTTY() const;                                            ///< Is the stream handled a TTY-like stream?

    StreamHandler& output();        ///< Output stream to use for all logging operations
    void setOutput(std::ostream*);  ///< Set the output stream

  private:
    explicit Logger(StreamHandler);  ///< Initialise a logging object

    std::vector<std::regex> allowed_exc_;  ///< List of enabled logging modules
    bool extended_{false};                 ///< Also print extra attributes?
    Level level_{Level::information};      ///< Logging threshold for the output stream
    StreamHandler output_{nullptr};        ///< Output stream to use for all logging operations
  };
}  // namespace cepgen::utils
namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const utils::Logger::Level&);
}

#define CG_LOG_MATCH(str, type) cepgen::utils::Logger::get().passExceptionRule(str, cepgen::utils::Logger::Level::type)
#define CG_LOG_LEVEL(type) cepgen::utils::Logger::get().setLevel(cepgen::utils::Logger::Level::type)

#endif
