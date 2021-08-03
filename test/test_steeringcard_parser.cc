#include <string>

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;

int main(int argc, char* argv[]) {
  string card;
  bool debug;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("card,i", "input card", &card)
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .parse();

  if (debug)
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debug;

  cepgen::initialise();

  try {
    CG_INFO("main") << "Parsing configuration from '" << card << ".";
    const auto* params = cepgen::card::Handler::parse(card);
    CG_INFO("main") << "Configuration parsed from '" << card << "':\n" << params;
  } catch (const cepgen::Exception& e) {
    e.dump();
    return -1;
  }
  return 0;
}
