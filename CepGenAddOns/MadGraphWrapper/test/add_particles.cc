#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();
  cepgen::initialise();

  const cepgen::pdgid_t my_part = 13;

  cepgen::PDG::get().define(cepgen::ParticleProperties(my_part, "la", "laurentino", 0., 42., 0., 3, true));
  CG_LOG << cepgen::PDG::get()(my_part);

  auto mg5 = cepgen::ProcessFactory::get().build(
      "mg5_aMC",
      cepgen::ParametersList()
          .set("kinematicsGenerator", cepgen::ParametersList().setName("coll2to4"s))
          .set("extraParticles", cepgen::ParametersList().set("la", cepgen::PDG::get()(my_part)))
          .set("process", "a a > la+ la-"s));

  return 0;
}
