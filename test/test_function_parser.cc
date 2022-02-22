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
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<string> parsers;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("parsers,p", "list of parsers to use", &parsers, vector<string>{})
      .parse();

  const double epsilon = 1.e-9;  // tolerance
  cepgen::initialise();

  if (parsers.empty())
    parsers = cepgen::utils::FunctionalFactory::get().modules();
  CG_LOG << "Will test with " << cepgen::utils::s("module", parsers.size(), true) << ": " << parsers;

  for (const auto& func : parsers) {
    CG_LOG << "Testing with \"" << func << "\" functional parser.";
    {  // test with a 1-variable function
      const double exp_result_test1 = 6.795704571;
      auto params = cepgen::ParametersList()
                        .set<std::string>("expression", "2.5*exp(0.1*x)")
                        .set<std::vector<std::string> >("variables", {"x"});
      auto test = cepgen::utils::FunctionalFactory::get().build(func, params);
      try {
        if (fabs((*test)(10.) - exp_result_test1) > epsilon) {
          CG_LOG << "Test 1.1 failed.";
          return -1;
        }
        if (fabs((*test)({10.}) - exp_result_test1) > epsilon) {
          CG_LOG << "Test 1.2 failed.";
          return -1;
        }
      } catch (const cepgen::Exception& e) {
        CG_LOG << "Test 1 failed.";
        return -1;
      }
      CG_LOG << "Test 1 passed.";
    }
    {  // test with an invalid function
      auto params = cepgen::ParametersList()
                        .set<std::string>("expression", "sqrt(x+x**3-log(10)")
                        .set<std::vector<std::string> >("variables", {"x"});
      try {
        auto test = cepgen::utils::FunctionalFactory::get().build(func, params);
        (*test)(10);
        CG_LOG << "Test 2 failed.";
        return -1;
      } catch (const cepgen::Exception& e) {
        CG_LOG << "Test 2 passed.";
      }
    }
    {  // test with a 2-variables function
      auto params = cepgen::ParametersList()
                        .set<std::string>("expression", "sqrt(a^2+b^2)")
                        .set<std::vector<std::string> >("variables", {"a", "b"});
      try {
        auto test = cepgen::utils::FunctionalFactory::get().build(func, params);
        if (fabs((*test)({3, 4}) - 5.0) > epsilon) {
          throw CG_ERROR("main") << "Exceeded the numerical limit.";
          return -1;
        }
        CG_LOG << "Test 3 passed.";
      } catch (const cepgen::Exception&) {
        CG_LOG << "Test 3 failed.";
        return -1;
      }
    }
    {  // test with an invalid function
      auto params = cepgen::ParametersList()
                        .set<std::string>("expression", "a***2")
                        .set<std::vector<std::string> >("variables", {"a"});
      try {
        auto test = cepgen::utils::FunctionalFactory::get().build(func, params);
        (*test)(10);
        (*test)({10});
        CG_LOG << "Test 4 failed";
        return -1;
      } catch (const cepgen::Exception& e) {
        CG_LOG << "Test 4 passed.";
      }
    }
  }

  return 0;
}
