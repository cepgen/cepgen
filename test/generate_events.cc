/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_card;
  int num_events;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("config,i", "path to the configuration file", &input_card, "Cards/lpair_cfg.py")
      .addOptionalArgument("num-events,n", "number of events to generate", &num_events, 10)
      .parse();

  cepgen::Generator gen;
  gen.setParameters(cepgen::card::Handler::parse(input_card));
  for (auto iev = 0; iev < num_events; ++iev)
    gen.next().dump();

  return 0;
}
