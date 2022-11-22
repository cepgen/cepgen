/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include <fstream>

#include "CepGen/Cards/Handler.h"
#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/ProgressBar.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Timer.h"

using namespace std;

int main(int argc, char* argv[]) {
  double num_sigma;
  string cfg_filename;
  string integrator;
  bool quiet;

  cepgen::ArgumentsParser argparse(argc, argv);
  argparse.addArgument("cfg,f", "configuration file", &cfg_filename)
      .addOptionalArgument("quiet,q", "quiet mode", &quiet, false)
      .addOptionalArgument("num-sigma,n", "max. number of std.dev.", &num_sigma, 3.)
      .addOptionalArgument("integrator,i", "type of integrator used", &integrator, "Vegas")
      .parse();

  if (quiet)
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::error;

  cepgen::utils::Timer tmr;
  cepgen::Generator gen;

  CG_LOG << "Testing with " << integrator << " integrator.";

  vector<string> failed_tests, passed_tests;

  CG_LOG << "Initial configuration time: " << tmr.elapsed() * 1.e3 << " ms.";
  tmr.reset();

  new cepgen::utils::AbortHandler;

  struct Test {
    string filename;
    double ref_cs, err_ref_cs;
  };
  vector<Test> tests;

  {  // parse the various tests to be performed
    ifstream cfg(cfg_filename);
    string line;
    while (!cfg.eof()) {
      getline(cfg, line);
      // skip the commented out lines
      line = cepgen::utils::ltrim(line);
      if (line.empty() || line[0] == '#')
        continue;
      stringstream os(line);
      Test test;
      os >> test.filename >> test.ref_cs >> test.err_ref_cs;
      tests.emplace_back(test);
      CG_DEBUG("main") << "Added test \"" << line << "\".";
    }
  }

  CG_LOG << "Will run " << cepgen::utils::s("test", tests.size()) << ".";

  std::unique_ptr<cepgen::utils::ProgressBar> progress;
  if (argparse.debugging())
    progress.reset(new cepgen::utils::ProgressBar(tests.size()));

  unsigned short num_tests = 0;
  for (const auto& test : tests) {
    const std::string filename = "test/physics/test_processes/" + test.filename + "_cfg.py";
    try {
      gen.parametersRef().clearProcess();
      gen.setParameters(cepgen::card::Handler::parse(filename));

      CG_DEBUG("main") << gen.parameters();

      gen.parametersRef().par_integrator.setName<std::string>(integrator);
      CG_LOG << "Process: " << gen.parameters()->processName() << "\n\t"
             << "File: " << filename << "\n\t"
             << "Configuration time: " << tmr.elapsed() * 1.e3 << " ms.";

      tmr.reset();

      double new_cs, err_new_cs;
      gen.computeXsection(new_cs, err_new_cs);

      const double ratio = new_cs / test.ref_cs;
      const double err_ratio = ratio * hypot(err_new_cs / new_cs, test.err_ref_cs / test.ref_cs);
      const double pull = (new_cs - test.ref_cs) / hypot(err_new_cs, test.err_ref_cs);

      const bool success = fabs(pull) < num_sigma;

      CG_LOG << "Computed cross section:\n\t"
             << "Ref.   = " << test.ref_cs << " +/- " << test.err_ref_cs << "\n\t"
             << "CepGen = " << new_cs << " +/- " << err_new_cs << "\n\t"
             << "Ratio: " << ratio << " +/- " << err_ratio << "\n\t"
             << "Pull: " << pull << " (abs(pull) " << (success ? "<" : ">") << " " << num_sigma << ").\n\t"
             << "Computation time: " << tmr.elapsed() * 1.e3 << " ms.";
      tmr.reset();

      const string test_res = cepgen::utils::format(
          "%-40s\tref=%g\tgot=%g\tratio=%g\tpull=%+10.5f", test.filename.c_str(), test.ref_cs, new_cs, ratio, pull);
      if (success)
        passed_tests.emplace_back(test_res);
      else
        failed_tests.emplace_back(test_res);
      ++num_tests;
      if (argparse.debugging())
        progress->update(num_tests);
      CG_LOG << "Test " << num_tests << "/" << tests.size() << " finished. "
             << "Success: " << cepgen::utils::yesno(success) << ".";
    } catch (const cepgen::Exception& e) {
      CG_LOG << "Test \"" << test.filename << "\" (located at " << filename << ") failed.";
      cepgen::Exception(e).dump();
    }
  }
  if (failed_tests.size() != 0) {
    ostringstream os_failed, os_passed;
    for (const auto& fail : failed_tests)
      os_failed << "  " << fail << endl;
    for (const auto& pass : passed_tests)
      os_passed << "  " << pass << endl;
    throw CG_FATAL("main") << "Some tests failed (abs(pull) > " << num_sigma << "):\n"
                           << os_failed.str() << "\n "
                           << "Passed tests:\n"
                           << os_passed.str() << ".";
  }

  CG_LOG << "ALL TESTS PASSED!";

  return 0;
}
