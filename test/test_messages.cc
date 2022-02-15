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

#include <fstream>

#include "CepGen/Utils/Logger.h"
#include "CepGen/Utils/Message.h"

using namespace cepgen;

int main() {
  //utils::Logger::get().level = utils::Logger::Level::debugInsideLoop;
  //utils::Logger::get().output = new std::ofstream("test.log");
  //utils::Logger::get().output = nullptr;

  //--- try with a bit of unicode too
  const std::string test_string = "Haha, ceci est un test à géométrie variable! ☺";
  for (int type = (int)LoggedMessage::MessageType::undefined; type <= (int)LoggedMessage::MessageType::warning;
       ++type) {
    LoggedMessage("Test", "", (LoggedMessage::MessageType)type) << test_string;
  }
  return 0;
}
