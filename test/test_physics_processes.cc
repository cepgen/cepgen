#include <cmath>
#include <fstream>
#include <iostream>

#include "CepGen/Cards/Handler.h"
#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/ProgressBar.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Timer.h"

using namespace std;
using namespace cepgen;

int main(int argc, char* argv[]) {
  double num_sigma;
  string cfg_filename;
  string integrator;
  bool debug, quiet;

  ArgumentsParser(argc, argv)
      .addArgument("cfg,c", "configuration file", &cfg_filename)
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .addOptionalArgument("quiet,q", "quiet mode", &quiet, false)
      .addOptionalArgument("num-sigma,n", "max. number of std.dev.", &num_sigma, 3.)
      .addOptionalArgument("integrator,i", "type of integrator used", &integrator, "Vegas")
      .parse();

  if (debug)
    utils::Logger::get().level = utils::Logger::Level::information;
  else if (quiet)
    utils::Logger::get().level = utils::Logger::Level::error;

  utils::Timer tmr;
  Generator gen;

  CG_LOG << "Testing with " << integrator << " integrator.";

  vector<string> failed_tests, passed_tests;

  CG_INFO("main") << "Initial configuration time: " << tmr.elapsed() * 1.e3 << " ms.";
  tmr.reset();

  new utils::AbortHandler;

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
      if (line[0] == '#' || line.empty())
        continue;
      stringstream os(line);
      Test test;
      os >> test.filename >> test.ref_cs >> test.err_ref_cs;
      tests.emplace_back(test);
    }
  }

  CG_INFO("main") << "Will run " << utils::s("test", tests.size()) << ".";

  std::unique_ptr<utils::ProgressBar> progress;
  if (debug)
    progress.reset(new utils::ProgressBar(tests.size()));

  try {
    unsigned short num_tests = 0;
    for (const auto& test : tests) {
      gen.parametersRef().clearProcess();

      const std::string filename = "test/test_processes/" + test.filename + "_cfg.py";
      gen.setParameters(cepgen::card::Handler::parse(filename));
      gen.parameters()->integrator->setName<std::string>(integrator);
      CG_INFO("main") << "Process: " << gen.parameters()->processName() << "\n\t"
                      << "File: " << filename << "\n\t"
                      << "Configuration time: " << tmr.elapsed() * 1.e3 << " ms.";

      tmr.reset();

      double new_cs, err_new_cs;
      gen.computeXsection(new_cs, err_new_cs);

      const double ratio = new_cs / test.ref_cs;
      const double err_ratio = ratio * hypot(err_new_cs / new_cs, test.err_ref_cs / test.ref_cs);
      const double pull = (new_cs - test.ref_cs) / hypot(err_new_cs, test.err_ref_cs);

      const bool success = fabs(pull) < num_sigma;

      CG_INFO("main") << "Computed cross section:\n\t"
                      << "Ref.   = " << test.ref_cs << " +/- " << test.err_ref_cs << "\n\t"
                      << "CepGen = " << new_cs << " +/- " << err_new_cs << "\n\t"
                      << "Ratio: " << ratio << " +/- " << err_ratio << "\n\t"
                      << "Pull: " << pull << " (abs(pull) " << (success ? "<" : ">") << " " << num_sigma << ").";

      CG_INFO("main") << "Computation time: " << tmr.elapsed() * 1.e3 << " ms.";
      tmr.reset();

      const string test_res = utils::format(
          "%-40s\tref=%g\tgot=%g\tratio=%g\tpull=%+10.5f", test.filename.c_str(), test.ref_cs, new_cs, ratio, pull);
      if (success)
        passed_tests.emplace_back(test_res);
      else
        failed_tests.emplace_back(test_res);
      ++num_tests;
      if (debug)
        progress->update(num_tests);
      CG_LOG << "Test " << num_tests << "/" << tests.size() << " finished. "
             << "Success: " << utils::yesno(success) << ".";
    }
  } catch (const Exception& e) {
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

  CG_INFO("main") << "ALL TESTS PASSED!";

  return 0;
}
