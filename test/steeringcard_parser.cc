/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021  Laurent Forthomme
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

#include <string>

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;

int main(int argc, char* argv[]) {
  string card;

  cepgen::initialise();

  cepgen::ArgumentsParser(argc, argv).addOptionalArgument("card,i", "input card", &card, "Cards/lpair_cfg.py").parse();

  try {
    CG_LOG << "Parsing configuration from '" << card << ".";
    const auto* params = cepgen::card::Handler::parseFile(card);
    CG_LOG << "Configuration parsed from '" << card << "':\n" << params;
  } catch (const cepgen::Exception& e) {
    e.dump();
    return -1;
  }
  return 0;
}
