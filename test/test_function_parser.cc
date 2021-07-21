#include <cmath>
#include <string>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"

using namespace std;

int main() {
  const double epsilon = 1.e-9;  // tolerance
  for (const auto& func : cepgen::utils::FunctionalFactory::get().modules()) {
    CG_LOG << "Testing with \"" << func << "\" functional parser.";
    {  // test with a 1-variable function
      const double exp_result_test1 = 6.795704571;
      auto test =
          cepgen::utils::FunctionalFactory::get().build(func,
                                                        cepgen::ParametersList()
                                                            .set<std::string>("expression", "2.5*exp(0.1*x)")
                                                            .set<std::vector<std::string> >("variables", {"x"}));
      if (fabs((*test)(10.) - exp_result_test1) > epsilon)
        throw CG_FATAL("main") << "Test 1.1 failed!";
      if (fabs((*test)({10.}) - exp_result_test1) > epsilon)
        throw CG_FATAL("main") << "Test 1.2 failed!";
      CG_LOG << "Test 1 passed!";
    }
    {  // test with an invalid function
      bool passed = false;
      try {
        auto test =
            cepgen::utils::FunctionalFactory::get().build(func,
                                                          cepgen::ParametersList()
                                                              .set<std::string>("expression", "sqrt(x+x**3-log(10)")
                                                              .set<std::vector<std::string> >("variables", {"x"}));
        (*test)(10);
      } catch (const cepgen::Exception& e) {
        CG_LOG << "Test 2 passed!";
        passed = true;
      }
      if (!passed)
        throw CG_FATAL("main") << "Test 2 failed!";
    }
    {  // test with a 2-variables function
      auto test =
          cepgen::utils::FunctionalFactory::get().build(func,
                                                        cepgen::ParametersList()
                                                            .set<std::string>("expression", "sqrt(a^2+b^2)")
                                                            .set<std::vector<std::string> >("variables", {"a", "b"}));
      if (fabs((*test)({3, 4}) - 5.0) > epsilon)
        throw CG_FATAL("main") << "Test 3 failed!";
      CG_LOG << "Test 3 passed!";
    }
    {  // test with an invalid function
      bool passed = true;
      try {
        auto test =
            cepgen::utils::FunctionalFactory::get().build(func,
                                                          cepgen::ParametersList()
                                                              .set<std::string>("expression", "a***2")
                                                              .set<std::vector<std::string> >("variables", {"a"}));
        (*test)(10);
        (*test)({10});
        passed = false;
      } catch (const cepgen::Exception& e) {
      }
      if (!passed)
        throw CG_FATAL("main") << "Test 4 failed!";
      CG_LOG << "Test 4 passed!";
    }
  }

  return 0;
}
