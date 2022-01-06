#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Modules/NamedModule.h"

int main() {
  class TestModule : public cepgen::NamedModule<std::string> {
  public:
    explicit TestModule(const cepgen::ParametersList& params) : cepgen::NamedModule<std::string>(params) {}
    static cepgen::ParametersDescription description() {
      cepgen::ParametersDescription desc("test_module");
      desc.add<int>("foo", 42);
      // description of a sub-collection of parameters
      cepgen::ParametersDescription submod("test_submodule");
      submod.add<double>("bar", 42.42);
      submod.add<std::string>("bat", "man").setDescription("What is in a 'bat'?");
      desc.add<cepgen::ParametersDescription>("sub_module_params", submod)
          .setDescription("A sub-collection of parameters");

      //desc.add<cepgen::ParametersList>("prout", cepgen::ParametersList());
      desc.add<std::string>("baz", "fourty-two").setDescription("A beautiful 'baz' name");
      return desc;
    }
  };

  CG_LOG << "Description of the test module:\n\n"
         << TestModule::description().describe()
         << "\nEquivalent parameters list: " << TestModule::description().parameters();

  return 0;
}
