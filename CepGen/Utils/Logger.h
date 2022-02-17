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
      /// Logging threshold for the output stream
      enum class Level { nothing = 0, error, warning, information, debug, debugInsideLoop };

    private:
      /// Initialize a logging object
      explicit Logger(std::ostream* os);

      std::vector<std::regex> allowed_exc_;
      bool extended_{false};

    public:
      /// Retrieve the running instance of the logger
      static Logger& get(std::ostream* os = nullptr);
      /// \brief Add a new rule to display exceptions/messages
      /// \param[in] rule Regex rule to handle
      void addExceptionRule(const std::string& rule);
      /// Collection of logging exceptions
      const std::vector<std::regex>& exceptionRules() const { return allowed_exc_; }
      /// Is the module set to be displayed/logged?
      /// \param[in] tmpl Module name to probe
      /// \param[in] lev Upper verbosity level
      bool passExceptionRule(const std::string& tmpl, const Level& lev) const;
      /// Also show extended information?
      bool extended() const { return extended_; }
      /// Set the extended information flag
      void setExtended(bool ext = true) { extended_ = ext; }

      /// Redirect the logger to a given output stream
      /// Logging threshold for the output stream
      Level level;
      /// Output stream to use for all logging operations
      std::ostream* output;
    };
  }  // namespace utils
  std::ostream& operator<<(std::ostream& os, const utils::Logger::Level&);
}  // namespace cepgen

#define CG_LOG_MATCH(str, type) cepgen::utils::Logger::get().passExceptionRule(str, cepgen::utils::Logger::Level::type)

#endif
