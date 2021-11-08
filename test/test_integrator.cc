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
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

using namespace std;

/// Generic process to test the integrator instance
class TestProcess : public cepgen::proc::Process {
public:
  /// Test process constructor
  inline explicit TestProcess(const string& func_mod, const string& formula, const vector<string>& args)
      : cepgen::proc::Process(cepgen::ParametersList().setName<string>(formula).set<string>("description", formula),
                              false),
        variables_(args.size()) {
    function_ = cepgen::utils::FunctionalFactory::get().build(
        func_mod, cepgen::ParametersList().set<string>("expression", formula).set<vector<string> >("arguments", args));
  }
  /// Process cloning method
  cepgen::proc::ProcessPtr clone() const override { return cepgen::proc::ProcessPtr(new TestProcess(*this)); }
  /// Dummy function to be called on phase space definition
  void prepareKinematics() override {
    for (auto& var : variables_)
      defineVariable(var, cepgen::proc::Process::Mapping::linear);
  }
  /// Dummy function to be called on events generation
  void fillKinematics(bool) override {}
  /// Generic formula to compute a weight out of a point in the phase space
  double computeWeight() override { return function_->operator()(variables_); }

private:
  vector<double> variables_;
  shared_ptr<cepgen::utils::Functional> function_;
};

int main(int argc, char* argv[]) {
  bool debug, quiet, run_all;
  double num_sigma;
  vector<string> integrators;
  string func_mod;

  cepgen::initialise();

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("num-sigma,n", "max. number of std.dev.", &num_sigma, 5.)
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .addOptionalArgument("integrator,i", "type of integrator used", &integrators, vector<string>{"Vegas"})
      .addOptionalArgument("functional,f", "type of functional parser user", &func_mod, "ROOT")
      .addOptionalArgument("all,a", "run the tests for all integrators", &run_all, false)
      .addOptionalArgument("quiet,q", "quiet mode", &quiet, false)
      .parse();

  if (debug)
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debug;
  else if (quiet)
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::nothing;
  else
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::information;

  //--- tests definition
  struct test_t {
    TestProcess process;
    double result;
    bool success;
  };

  vector<test_t> tests = {
      {TestProcess(func_mod, "x^2+y^2", {"x", "y"}), 2. / 3, false},
      {TestProcess(func_mod, "x+y^2+z^3", {"x", "y", "z"}), 13. / 12., false},
      {TestProcess(func_mod, "1./(1.-cos(x*3.141592654)*cos(y*3.141592654)*cos(z*3.141592654))", {"x", "y", "z"}),
       1.3932039296856768591842462603255,
       false},
  };

  //--- integrator definition
  if (run_all)
    // will perform the test with all integrators
    integrators = cepgen::IntegratorFactory::get().modules();

  CG_LOG << "Will test with " << cepgen::utils::s("integrator", integrators.size(), true) << ": " << integrators;

  cepgen::Parameters params;

  for (const auto& integrator : integrators) {
    CG_LOG << "Running with " << integrator << " integrator.";
    auto integr = cepgen::IntegratorFactory::get().build(integrator);

    //--- integration part
    size_t i = 0;
    double result, error;
    for (auto& test : tests) {
      params.setProcess(test.process.clone());
      cepgen::ProcessIntegrand integrand(&params);
      integr->setIntegrand(integrand);
      integr->integrate(result, error);
      test.success = error / result < 1.e-6 || (fabs(test.result - result) <= num_sigma * error);
      if (debug)
        CG_LOG << "Test " << i << ": ref.: " << test.result << ", result: " << result << " +/- " << error << ".";
      ++i;
    }

    bool success = true;
    i = 0;
    for (const auto& test : tests) {
      CG_LOG << "Test " << i++ << " passed: " << cepgen::utils::yesno(test.success);
      success &= test.success;
    }
    if (!success)
      throw CG_FATAL("main") << integrator << " integrator tests failed!";
  }
  return 0;
}
