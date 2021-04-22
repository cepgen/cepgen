#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Parameters.h"

#include "CepGen/Utils/ArgumentsParser.h"

#include <string>

using namespace std;
using namespace cepgen;

int main(int argc, char* argv[]) {
  string card;

  ArgumentsParser(argc, argv).addArgument("card,i", "input card", &card).parse();

  try {
    CG_INFO("main") << card::PythonHandler(ParametersList().set<string>("filename", card)).parameters();
  } catch (const Exception& e) {
    e.dump();
  }
  return 0;
}
