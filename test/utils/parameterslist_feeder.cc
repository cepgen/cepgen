#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  {
    const string fed = "test/of/key:value";
    cepgen::ParametersList plist;
    plist.feed(fed);
    CG_TEST_EQUAL(plist.get<cepgen::ParametersList>("test").get<cepgen::ParametersList>("of").get<std::string>("key"),
                  "value",
                  "parameters list chain");
    CG_DEBUG("main") << "Resulting parameters list: " << plist << ".";
    CG_TEST_EQUAL(plist.serialise(), fed, "parameters list serialisation");
  }
  {
    cepgen::ParametersList plist;
    plist.feed("foo:3.14").feed("bar:2").feed("baz:2e3");
    CG_DEBUG("main") << "Resulting parameters list: " << plist << ".";
    CG_TEST_EQUAL(plist.get<double>("foo"), 3.14, "float parsing");
    CG_TEST_EQUAL(plist.get<int>("bar"), 2, "integer parsing");
    CG_TEST_EQUAL(plist.get<double>("baz"), 2000., "float (from Ee notation");
    plist.feed("bat:5E10").feed("foo:42");
    CG_TEST_EQUAL(plist.get<double>("bat"), 5e10, "float (from re-parsing)");
    CG_TEST_EQUAL(plist.get<int>("foo"), 42, "integer (from re-parsing)");
    CG_TEST_EQUAL(plist.get<double>("foo", -1.), -1., "integer as float (from re-parsing)");
  }
  {
    auto fed = "this/is/a:test,this/works:true,that/{one:42,other:3.141592}";
    cepgen::ParametersList plist;
    plist.feed(fed);
    const auto re_fed = cepgen::ParametersList().feed(plist.serialise());
    CG_DEBUG("main") << "\n"
                     << "Fed string: " << fed << "\n"
                     << "Fed parameters list:\n"
                     << cepgen::ParametersDescription(plist) << "\n"
                     << "Re-serialised string: " << plist.serialise() << "\n"
                     << "Re-fed parameters list:\n"
                     << cepgen::ParametersDescription(re_fed) << ".\n"
                     << "Diff:\n"
                     << cepgen::ParametersDescription(plist.diff(re_fed)) << ".";
    CG_TEST_EQUAL(re_fed, plist, "serialised parameters list parsing");
  }
  CG_TEST_EXCEPT([]() { cepgen::ParametersList().feed("invalid/string/{{fed:true}"); }, "parsing of an invalid string");

  CG_TEST_SUMMARY;
}
