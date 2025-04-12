#include <Pythia8/Pythia.h>

#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main() {
  cepgen::initialise();

  constexpr auto seed1 = 1234567, seed2 = 7654321;

  const auto cg_pythia = cepgen::EventModifierFactory::get().build("pythia8");
  cg_pythia->readString("Random:seed = " + to_string(seed1));

  const auto pythia = cg_pythia->engine<Pythia8::Pythia>();
  CG_TEST_EQUAL(pythia->checkVersion(), true, "Pythia 8 object version");
  CG_TEST_EQUAL(pythia->mode("Random:seed"), seed1, "Parameter set on wrapper");

  pythia->readString("Random:seed = " + to_string(seed2));
  CG_TEST_EQUAL(pythia->mode("Random:seed"), seed2, "Parameter set on engine");

  CG_TEST_SUMMARY;
}
