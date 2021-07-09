#include <string>

#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;
using namespace cepgen;

int main(int argc, char* argv[]) {
  string card;
  bool debug;

  ArgumentsParser(argc, argv)
      .addOptionalArgument("card,i", "input card", &card, "Cards/lpair_cfg.py")
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .parse();

  if (debug)
    utils::Logger::get().level = utils::Logger::Level::debug;

  try {
    CG_INFO("main") << "Parsing configuration from '" << card << ".";
    const card::PythonHandler parsed(ParametersList().set<string>("filename", card));
    CG_INFO("main") << "Configuration parsed from '" << card << "':\n" << parsed.runtimeParameters();
  } catch (const Exception& e) {
    e.dump();
  }
  return 0;
}
