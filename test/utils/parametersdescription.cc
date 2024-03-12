#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  class TestModule : public cepgen::SteeredObject<TestModule> {
  public:
    explicit TestModule(const cepgen::ParametersList& params) : SteeredObject(params) {}
    static cepgen::ParametersDescription description() {
      cepgen::ParametersDescription desc("test_module");
      desc.add<int>("foo", 42);
      {
        cepgen::ParametersDescription submod("test_submodule");
        submod.add<double>("bar", 42.42);
        submod.add<std::string>("bat", "man").setDescription("What is in a 'bat'?");
        desc.add("sub_module_params", submod).setDescription("A sub-collection of parameters");
      }
      desc.add<std::string>("baz", "fourty-two").setDescription("A beautiful 'baz' name");
      return desc;
    }
  };

  const auto mod = TestModule(cepgen::ParametersList().set<int>("foo", 21));
  CG_DEBUG("main") << "Description of the test module:\n"
                   << TestModule::description().describe()
                   << "\nEquivalent parameters list: " << TestModule::description().parameters()
                   << "\nSteered test module:\n"
                   << cepgen::ParametersDescription(mod.parameters());
  CG_TEST_EQUAL(mod.parameters().get<std::string>("baz"), "fourty-two", "un-steered parameter in module");
  CG_TEST_EQUAL(mod.parameters().get<int>("foo"), 21, "steered parameter in module");
  CG_TEST_EQUAL(mod.parameters().get<cepgen::ParametersList>("sub_module_params").get<std::string>("bat"),
                "man",
                "un-steered parameters in module's sub-parameters");

  CG_TEST_SUMMARY;
}
