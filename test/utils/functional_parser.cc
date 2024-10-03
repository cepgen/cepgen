/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <string>

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::initialise();

  vector<string> parsers;
  bool verbose;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("parsers,p", "list of parsers to use", &parsers, cepgen::FunctionalFactory::get().modules())
      .addOptionalArgument("verbose", "verbose mode", &verbose, false)
      .parse();

  CG_TEST_DEBUG(verbose);
  constexpr double epsilon = 1.e-9;  // tolerance

  CG_LOG << "Will test with " << cepgen::utils::s("module", parsers.size(), true) << ": " << parsers;

  for (const auto& func : parsers) {
    CG_LOG << "Testing with \"" << func << "\" functional parser.";
    {  // test with a 1-variable function
      constexpr double exp_result_test1 = 6.795704571;
      CG_LOG << cepgen::utils::Functional::fromExpression("2.5*exp(0.1*x)", {"x"});
      auto test = cepgen::FunctionalFactory::get().build(
          func, cepgen::utils::Functional::fromExpression("2.5*exp(0.1*x)", {"x"}));
      CG_TEST(fabs((*test)(10.) - exp_result_test1) <= epsilon, "single argument functional");
      CG_TEST(fabs((*test)({10.}) - exp_result_test1) <= epsilon, "multiple-argument functional");
    }
    {  // test with an invalid function
      auto test_invalid = [&func]() {
        auto test = cepgen::FunctionalFactory::get().build(
            func, cepgen::utils::Functional::fromExpression("sqrt(x+x**3-log(10)", {"x"}));
        (void)(*test)(10);
      };
      CG_TEST_EXCEPT(test_invalid, "invalid function parsing");
    }
    {  // test with a 2-variables function
      try {
        auto test = cepgen::FunctionalFactory::get().build(
            func, cepgen::utils::Functional::fromExpression("sqrt(a^2+b^2)", {"a", "b"}));
        CG_TEST(fabs((*test)({3, 4}) - 5.0) <= epsilon, "two-variables function");
      } catch (const cepgen::Exception&) {
        CG_LOG << "Test 3 failed.";
        return -1;
      }
    }
    {  // test with an invalid function
      try {
        auto test =
            cepgen::FunctionalFactory::get().build(func, cepgen::utils::Functional::fromExpression("a***2", {"a"}));
        (void)(*test)(10);
        (void)(*test)({10});
        CG_LOG << "Test 4 failed";
        return -1;
      } catch (const cepgen::Exception&) {
        CG_LOG << "Test 4 passed.";
      }
    }
  }
  CG_TEST_SUMMARY;
}
