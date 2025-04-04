/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2015-2025  Laurent Forthomme
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

#include <iostream>

#include "CepGen/Utils/Logger.h"
#include "CepGen/Utils/Message.h"

using namespace cepgen::utils;

Logger::Logger(StreamHandler os) : output_(std::move(os)) {
  if (output_.get() != &std::cout)
    CG_INFO("Logger") << "New logger initialised for output@0x" << output_.get() << ".";
}

Logger& Logger::get(std::ostream* os) {
  static Logger log(os ? std::unique_ptr<std::ostream>(os) : StreamHandler(&std::cout, [](std::ostream*) -> void {}));
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

void Logger::setOutput(std::ostream* os) { output_.reset(os); }

Logger::StreamHandler& Logger::output() {
  static StreamHandler empty_output{nullptr};
  if (level_ == Level::nothing)
    return empty_output;
  return output_;
}

bool Logger::isTTY() const { return output_.get() == &std::cout || output_.get() == &std::cerr; }

namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const Logger::Level& level) {
    switch (level) {
      case Logger::Level::nothing:
        return os << "None";
      case Logger::Level::error:
        return os << "Errors";
      case Logger::Level::warning:
        return os << "Warnings";
      case Logger::Level::information:
        return os << "Infos";
      case Logger::Level::debug:
        return os << "Debug";
      case Logger::Level::debugInsideLoop:
        return os << "Debug (in loops)";
    }
    return os << "Unknown severity";
  }
}  // namespace cepgen
