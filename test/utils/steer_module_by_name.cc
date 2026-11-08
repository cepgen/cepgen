#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();
  {
    const string fed = "text<test/of/key:value<width:40";
    auto drawer = cepgen::DrawerFactory::get().build(fed);
    CG_TEST_EQUAL(drawer->parameters().get<bool>("colourise"), true, "unaffected parameter");
    CG_TEST_EQUAL(drawer->parameters().get<cepgen::ParametersList>("test"),
                  cepgen::ParametersList().feed("of/key:value"),
                  "unrelated parameters list hierarchy");
    CG_TEST_EQUAL(drawer->parameters().get<int>("width"), 40, "extra integer parameter");
  }

  CG_TEST_SUMMARY;
}
