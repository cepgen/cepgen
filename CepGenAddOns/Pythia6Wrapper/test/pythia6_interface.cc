#include "CepGen/Generator.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenAddOns/Pythia6Wrapper/Pythia6Interface.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::initialise();

  for (auto pdgid : {2212, 11, 13, 22}) {
    ostringstream os;
    os << cepgen::PDG::get()(pdgid).name << "/" << pythia6::pyname(pdgid);
    CG_TEST_EQUIV(pythia6::pymass(pdgid), cepgen::PDG::get().mass(pdgid), os.str() + " mass");
  }

  CG_TEST_SUMMARY;
}
