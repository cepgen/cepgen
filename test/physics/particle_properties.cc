#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::ParametersList plist;
  plist.set<string>("name", "laurenteron");

  cepgen::ParticleProperties prop(plist);
  CG_TEST_EQUAL(prop.name, plist.get<string>("name"), "custom particle name");

  prop.pdgid = 42;
  CG_TEST_EQUAL(prop.parameters().get<cepgen::pdgid_t>("pdgid"), prop.pdgid, "post-defined particle id change");

  CG_TEST_SUMMARY;
}
