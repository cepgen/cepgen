/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  CG_DEBUG("main") << "Will build a Python cards handler.";

  auto card = cepgen::CardsHandlerFactory::get().build(".py");
  card->parseString(R"(
from Config.PDG_cfi import PDG, registerParticle
registerParticle(name='teston', pdgid=42, mass=42.42, width=1.1))",
                    nullptr);

  CG_DEBUG("main") << "Configuration string successfully parsed.";

  const auto& teston = cepgen::PDG::get()(42);
  CG_TEST_EQUAL(teston.pdgid, 42, "new particle PDG id");
  CG_TEST_EQUAL(teston.name, "teston", "new particle name");
  CG_TEST_EQUAL(teston.mass, 42.42, "new particle mass");
  CG_TEST_EQUAL(teston.width, 1.1, "new particle width");

  CG_TEST_SUMMARY;
}
