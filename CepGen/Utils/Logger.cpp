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

#include <iostream>  // for cout

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Logger.h"

namespace cepgen {
  namespace utils {
    Logger::Logger(std::ostream* os) : output_(os) {
      if (output_ != &std::cout)
        CG_INFO("Logger") << "New logger initialised for output@" << output_ << ".";
    }

    Logger& Logger::get(std::ostream* os) {
      static Logger log(os ? os : &std::cout);
      return log;
    }

    void Logger::addExceptionRule(const std::string& rule) {
      allowed_exc_.emplace_back(rule, std::regex_constants::extended);
    }

    bool Logger::passExceptionRule(const std::string& tmpl, const Level& lev) const {
      if (level_ >= lev)
        return true;
      if (allowed_exc_.empty())
        return false;
      for (const auto& rule : allowed_exc_)
        try {
          if (std::regex_match(tmpl, rule))
            return true;
        } catch (const std::regex_error& err) {
          throw std::runtime_error("Failed to evaluate regex for logging tool.\n" + std::string(err.what()));
        }
      return false;
    }

    std::ostream* Logger::output() {
      static std::ostream empty_output(nullptr);
      if (level_ == Level::nothing)
        return &empty_output;
      return output_;
    }
  }  // namespace utils

  std::ostream& operator<<(std::ostream& os, const utils::Logger::Level& lvl) {
    switch (lvl) {
      case utils::Logger::Level::nothing:
        return os << "None";
      case utils::Logger::Level::error:
        return os << "Errors";
      case utils::Logger::Level::warning:
        return os << "Warnings";
      case utils::Logger::Level::information:
        return os << "Infos";
      case utils::Logger::Level::debug:
        return os << "Debug";
      case utils::Logger::Level::debugInsideLoop:
        return os << "Debug (in loops)";
    }
    return os;
  }
}  // namespace cepgen
