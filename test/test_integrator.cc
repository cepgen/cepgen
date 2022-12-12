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

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Integration/FunctionalIntegrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  bool quiet;
  double num_sigma;
  vector<string> integrators;
  string func_mod;

  cepgen::initialise();
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("num-sigma,n", "max. number of std.dev.", &num_sigma, 5.)
      .addOptionalArgument("integrator,i",
                           "type of integrator used",
                           &integrators,
                           cepgen::IntegratorFactory::get().modules())  // by default, all integrators are tested
      .addOptionalArgument("functional,f", "type of functional parser user", &func_mod, "ROOT")
      .addOptionalArgument("quiet,q", "quiet mode", &quiet, false)
      .parse();

  if (quiet)
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::nothing;

  //--- tests definition
  struct test_t {
    cepgen::FunctionalIntegrand integrand;
    double result;
    bool success{false};
  };

  vector<test_t> tests;
  tests.emplace_back(test_t{cepgen::FunctionalIntegrand("x^2+y^2", {"x", "y"}, func_mod), 2. / 3});
  tests.emplace_back(test_t{cepgen::FunctionalIntegrand("x+y^2+z^3", {"x", "y", "z"}, func_mod), 13. / 12.});
  tests.emplace_back(
      test_t{cepgen::FunctionalIntegrand(
                 "1./(1.-cos(x*3.141592654)*cos(y*3.141592654)*cos(z*3.141592654))", {"x", "y", "z"}, func_mod),
             1.3932039296856768591842462603255});

  CG_LOG << "Will test with " << cepgen::utils::s("integrator", integrators.size(), true) << ": " << integrators;

  for (const auto& integrator : integrators) {
    auto integr = cepgen::IntegratorFactory::get().build(integrator);

    //--- integration part
    size_t i = 0;
    double result, error;
    for (auto& test : tests) {
      integr->setIntegrand(test.integrand);
      integr->integrate(result, error);
      CG_DEBUG("main") << "Test " << i << ": ref.: " << test.result << ", result: " << result << " +/- " << error
                       << ".";
      CG_TEST(error / result < 1.e-6 || (fabs(test.result - result) <= num_sigma * error),
              integrator + " test " + std::to_string(i));
      ++i;
    }
  }
  CG_TEST_SUMMARY;
}
