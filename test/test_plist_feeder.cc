#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  {
    const string feeded = "test/of/key=value";
    cepgen::ParametersList plist;
    plist.feed(feeded);
    if (plist.get<cepgen::ParametersList>("test").get<cepgen::ParametersList>("of").get<std::string>("key") !=
        "value") {
      CG_LOG << "Failed to parse a parameters lists chain. Result=" << plist << ".";
      return -1;
    }
    CG_DEBUG("main") << "Resulting parameters list: " << plist << ".";
    if (plist.serialise() != feeded) {
      CG_LOG << "Failed to serialise the parameters list. Result=" << plist.serialise() << ".";
      return -1;
    }
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
  {
    auto feeded = "this/is/a=test,this/works=true,that/{one=42,other=3.141592}";
    cepgen::ParametersList plist;
    plist.feed(feeded);
    CG_DEBUG("main") << "\n\t"
                     << "Feeded string: " << feeded << "\n\t"
                     << "Fed parameters list: " << cepgen::ParametersDescription(plist) << "\n\t"
                     << "Re-serialised string: " << plist.serialise();
    if (cepgen::ParametersList().feed(plist.serialise()) != plist) {
      CG_LOG << "Failed to parse a serialised parameters list. Result=" << plist.serialise() << ".";
      return -1;
    }
  }
  {  // test with an invalid string feeded
    auto feeded = "invalid/string/{{feeded=true}";
    cepgen::ParametersList plist;
    try {
      plist.feed(feeded);
      CG_LOG << "Failed test of parsing an invalid string: " << feeded << ".";
      return -1;
    } catch (const cepgen::Exception&) {
      CG_DEBUG("main") << "Passed invalid string feeder test.";
    }
  }

  return 0;
}
