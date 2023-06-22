/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2015-2023  Laurent Forthomme
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

namespace cepgen {
  namespace utils {
    /// \brief General purposes logger
    /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
    /// \date 15 Oct 2015
    class Logger {
    public:
      /// Retrieve the running instance of the logger
      static Logger& get(std::ostream* = nullptr);

      /// \brief Add a new rule to display exceptions/messages
      /// \param[in] rule Regex rule to handle
      void addExceptionRule(const std::string& rule);
      /// Collection of logging exceptions
      const std::vector<std::regex>& exceptionRules() const { return allowed_exc_; }

      typedef std::unique_ptr<std::ostream, std::function<void(std::ostream*)> > StreamHandler;

      /// Logging threshold for the output stream
      enum class Level { nothing = 0, error, warning, information, debug, debugInsideLoop };
      /// Is the module set to be displayed/logged?
      /// \param[in] tmpl Module name to probe
      /// \param[in] lev Upper verbosity level
      bool passExceptionRule(const std::string& tmpl, const Level& lev) const;
      /// Logging threshold
      Level level() const { return level_; }
      /// Set the logging threshold
      void setLevel(Level level) { level_ = level; }
      /// Also show extended information?
      bool extended() const { return extended_; }
      /// Set the extended information flag
      void setExtended(bool ext = true) { extended_ = ext; }
      /// Is the stream handled a TTY-like stream?
      bool isTTY() const;

      /// Output stream to use for all logging operations
      StreamHandler& output();
      /// Set the output stream
      void setOutput(std::ostream*);

    private:
      /// Initialize a logging object
      explicit Logger(StreamHandler);

      /// List of enabled logging modules
      std::vector<std::regex> allowed_exc_;
      /// Also print extra attributes?
      bool extended_{false};
      /// Logging threshold for the output stream
      Level level_{Level::information};
      /// Output stream to use for all logging operations
      StreamHandler output_{nullptr};
    };
  }  // namespace utils
  std::ostream& operator<<(std::ostream& os, const utils::Logger::Level&);
}  // namespace cepgen

#define CG_LOG_MATCH(str, type) cepgen::utils::Logger::get().passExceptionRule(str, cepgen::utils::Logger::Level::type)
#define CG_LOG_LEVEL(type) cepgen::utils::Logger::get().setLevel(cepgen::utils::Logger::Level::type)

#endif
