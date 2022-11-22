#include <cmath>

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"

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

  if (test.parameters().get<int>("foo") != 42 - 19) {
    CG_LOG << "Test with int failed";
    return -1;
  } else if (test.parameters().get<double>("bar") != 2. * M_PI) {
    CG_LOG << "Test with float failed";
    return -1;
  } else if (test.parameters().get<std::string>("baz") != "☺") {
    CG_LOG << "Test with string failed";
    return -1;
  } else if (!test.parameters().get<bool>("bat")) {
    CG_LOG << "Test with boolean failed";
    return -1;
  }

  test.setParameters(cepgen::ParametersList().set<int>("foo", 41));
  if (test.parameters().get<int>("foo") != 41) {
    CG_LOG << "Test with parameters object-set int failed";
    return -1;
  }

  test.foo = 45;
  if (test.parameters().get<int>("foo") != 45) {
    CG_LOG << "Test with object-set int failed";
    return -1;
  }

  CG_LOG << "All tests passed.";

  return 0;
}
