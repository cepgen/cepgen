#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::ParametersList plist;
  plist.set<string>("name", "laurenteron");

  cepgen::ParticleProperties prop(plist);
  if (prop.name != plist.get<string>("name")) {
    CG_LOG << "Failed to specify particle name: " << prop.name << ".";
    return -1;
  }

  prop.pdgid = 42;
  if (prop.parameters().get<cepgen::pdgid_t>("pdgid") != prop.pdgid) {
    CG_LOG << "Failed to retrieve particle id from plist once specified in object: " << prop.parameters() << ".";
    return -1;
  }

  CG_LOG << prop;
  CG_LOG << prop.parameters();

  return 0;
}
