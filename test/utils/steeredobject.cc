#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std::string_literals;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  struct TestContainer : cepgen::SteeredObject<TestContainer> {
    explicit TestContainer() {
      (*this).add("foo", foo).add("bar", bar).add("baz", baz).add("bat", bat).add("ban", ban);
    }
    static cepgen::ParametersDescription description() {
      auto desc = cepgen::ParametersDescription();
      desc.add("foo", 42);
      desc.add("bar", M_PI);
      desc.add("baz", "test™"s);
      desc.add("bat", false);
      return desc;
    }
    int foo;
    double bar;
    std::string baz;
    bool bat;
    cepgen::pdgid_t ban;
  };

  TestContainer test;

  CG_TEST_EQUAL(test.foo, 42, "integer retrieval from members");
  CG_TEST_EQUAL(test.bar, M_PI, "float retrieval from members");
  CG_TEST_EQUAL(test.baz, "test™", "string retrieval from members");
  CG_TEST_EQUAL(test.bat, false, "boolean retrieval from members");
  CG_TEST_EQUAL(test.ban, 0, "PDG id retrieval from members");

  test.foo -= 19;
  test.bat = true;
  test.bar *= 2.;
  test.baz = "☺";
  test.ban = cepgen::PDG::photon;

  CG_TEST_EQUAL(test.parameters().get<int>("foo"), 42 - 19, "integer retrieval from parameters");
  CG_TEST_EQUAL(test.parameters().get<double>("bar"), 2. * M_PI, "float retrieval from parameters");
  CG_TEST_EQUAL(test.parameters().get<std::string>("baz"), "☺", "string retrieval from parameters");
  CG_TEST_EQUAL(test.parameters().get<bool>("bat"), true, "boolean retrieval from parameters");
  CG_TEST_EQUAL(test.parameters().get<cepgen::pdgid_t>("ban"), cepgen::PDG::photon, "PDG id retrieval from parameters");

  test.setParameters(cepgen::ParametersList().set<int>("foo", 41));
  CG_TEST_EQUAL(test.foo, 41, "integer retrieval from parameters-set object");
  test.setParameters(cepgen::ParametersList().set<cepgen::pdgid_t>("ban", cepgen::PDG::gluon));
  CG_TEST_EQUAL(test.ban, cepgen::PDG::gluon, "PDG retrieval from parameters-set object");

  test.foo = 45;
  CG_TEST_EQUAL(test.parameters().get<int>("foo"), 45, "integer retrieval from object-set parameters");

  CG_TEST_SUMMARY;
}
