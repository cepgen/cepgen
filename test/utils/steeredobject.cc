#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  struct TestContainer : cepgen::SteeredObject<TestContainer> {
    explicit TestContainer() { (*this).add("foo", foo).add("bar", bar).add("baz", baz).add("bat", bat); }
    static cepgen::ParametersDescription description() {
      auto desc = cepgen::ParametersDescription();
      desc.add<int>("foo", 42);
      desc.add<double>("bar", M_PI);
      desc.add<std::string>("baz", __PRETTY_FUNCTION__);
      desc.add<bool>("bat", false);
      return desc;
    }
    int foo;
    double bar;
    std::string baz;
    bool bat;
  };

  TestContainer test;
  test.foo -= 19;
  test.bat = true;
  test.bar *= 2.;
  test.baz = "☺";

  CG_TEST(test.parameters().get<int>("foo") == 42 - 19, "integer retrieval from parameters");
  CG_TEST(test.parameters().get<double>("bar") == 2. * M_PI, "float retrieval from parameters");
  CG_TEST(test.parameters().get<std::string>("baz") == "☺", "string retrieval from parameters");
  CG_TEST(test.parameters().get<bool>("bat") == true, "boolean retrieval from parameters");

  test.setParameters(cepgen::ParametersList().set<int>("foo", 41));
  CG_TEST(test.foo == 41, "integer retrieval from parameters-set object");

  test.foo = 45;
  CG_TEST(test.parameters().get<int>("foo") == 45, "integer retrieval from object-set parameters");


  CG_TEST_SUMMARY;
}
