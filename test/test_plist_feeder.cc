#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  {
    cepgen::ParametersList plist;
    plist.feed("test/of/key=value");
    if (plist.get<cepgen::ParametersList>("test").get<cepgen::ParametersList>("of").get<std::string>("key") !=
        "value") {
      CG_LOG << "Failed to parse a parameters lists chain. Result=" << plist << ".";
      return -1;
    }
    CG_DEBUG("main") << "Resulting parameters list: " << plist << ".";
  }
  {
    cepgen::ParametersList plist;
    plist.feed("foo=3.14").feed("bar=2").feed("baz=2e3");
    if (plist.get<double>("foo") != 3.14 || plist.get<int>("bar") != 2 || plist.get<double>("baz") != 2000.) {
      CG_LOG << "Failed to parse an integer/float parameters list. Result=" << plist << ".";
      return -1;
    }
    plist.feed("bat=5E10").feed("foo=42");
    if (plist.get<double>("bat") != 5.e10 || plist.get<int>("foo") != 42) {
      CG_LOG << "Failed to re-parse an integer/float parameters list. Result=" << plist << ".";
      return -1;
    }
    CG_DEBUG("main") << "Resulting parameters list: " << plist << ".";
  }

  return 0;
}
