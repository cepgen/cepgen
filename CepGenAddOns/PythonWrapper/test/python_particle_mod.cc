#include "CepGen/Cards/Handler.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenAddOns/PythonWrapper/Environment.h"

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  auto card = cepgen::CardsHandlerFactory::get().build(".py");
  CG_DEBUG("main") << "Python cards handler successfully constructed.";

  card->parseString(R"(
from Config.PDG_cfi import PDG, registerParticle
registerParticle(name='teston', pdgid=42, mass=42.42, width=1.1))",
                    new cepgen::Parameters);
  CG_DEBUG("main") << "Configuration string successfully parsed.";

  const auto& teston = cepgen::PDG::get()(42);
  CG_TEST_EQUAL(teston.pdgid, 42, "new particle PDG id");
  CG_TEST_EQUAL(teston.name, "teston", "new particle name");
  CG_TEST_EQUAL(teston.mass, 42.42, "new particle mass");
  CG_TEST_EQUAL(teston.width, 1.1, "new particle width");

  CG_TEST_SUMMARY;
}
