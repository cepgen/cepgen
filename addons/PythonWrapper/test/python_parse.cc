/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  auto card = cepgen::CardsHandlerFactory::get().build(".py");
  card->parseCommands({R"(import Config.Core as cepgen
lim = (42.42, 420.420)
vec_lim = [(0., 1.), (1., 2.)]
)"});

  // retrieve all parameters parsed from Python commands
  const auto& parsed_params = card->parameters().get<cepgen::ParametersList>("parsed");

  const auto lim = cepgen::Limits{42.42, 420.420};
  CG_TEST_EQUAL(parsed_params.get<cepgen::Limits>("lim"), lim, "limits");
  const vector<cepgen::Limits> vec_lim = {{0., 1.}, {1., 2.}};
  CG_TEST_EQUAL(parsed_params.get<vector<cepgen::Limits> >("vec_lim"), vec_lim, "vector of limits");

  CG_TEST_SUMMARY;
}
