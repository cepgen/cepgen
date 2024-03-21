/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Test.h"
#include "CepGen/Utils/Timer.h"

using namespace std;

int main(int argc, char* argv[]) {
  double num_sigma;
  string cfg_filename;
  string integrator;
  bool verbose, quiet;

  auto args = cepgen::ArgumentsParser(argc, argv)
                  .addOptionalArgument("cfg,f",
                                       "configuration file",
                                       &cfg_filename,
                                       cepgen::utils::env::get("CEPGEN_PATH") + "/test/physics/test_processes.cfg")
                  .addOptionalArgument("verbose,v", "verbose mode", &verbose, false)
                  .addOptionalArgument("quiet,q", "quiet mode", &quiet, false)
                  .addOptionalArgument("num-sigma,n", "max. number of std.dev.", &num_sigma, 3.)
                  .addOptionalArgument("integrator,i", "type of integrator used", &integrator, "MISER")
                  .parse();

  if (!args.debugging() && !verbose)
    CG_LOG_LEVEL(warning);

  cepgen::utils::Timer tmr;
  cepgen::Generator gen;

  CG_TEST_DEBUG(verbose);
  if (quiet)
    CG_LOG_LEVEL(warning);

  new cepgen::utils::AbortHandler;

  struct Test {
    string filename;
    cepgen::Value ref_cs;
  };
  vector<Test> tests;

  {  // parse the various tests to be performed
    ifstream cfg(cfg_filename);
    string line;
    while (getline(cfg, line)) {
      // skip the commented out lines
      line = cepgen::utils::ltrim(line);
      if (line.empty() || line[0] == '#')
        continue;
      Test test;
      double ref, err;
      stringstream os(line);
      os >> test.filename >> ref >> err;
      test.ref_cs = cepgen::Value{ref, err};
      tests.emplace_back(test);
      CG_DEBUG("main") << "Added test '" << test.filename << "' with expected cross section: " << test.ref_cs << " pb.";
    }
  }

  CG_LOG << "Will run " << cepgen::utils::s("test", tests.size()) << " with " << integrator << " integrator.";
  CG_LOG << "Initial configuration time: " << tmr.elapsed() * 1.e3 << " ms.";
  tmr.reset();

  for (const auto& test : tests) {
    const std::string filename = "TestProcesses/" + test.filename + "_cfg.py";
    try {
      gen.parseRunParameters(filename);
      gen.runParameters().integrator() = cepgen::IntegratorFactory::get().describeParameters(integrator).parameters();
      CG_DEBUG("main") << "Process: " << gen.runParameters().processName() << "\n\t"
                       << "File: " << filename << "\n\t"
                       << "Configuration time: " << tmr.elapsed() * 1.e3 << " ms.";
      tmr.reset();

      const auto new_cs = gen.computeXsection();
      const auto ratio = new_cs / test.ref_cs,
                 pull = (new_cs - test.ref_cs) / hypot(new_cs.uncertainty(), test.ref_cs.uncertainty());

      CG_DEBUG("main") << "Computed cross section:\n\t"
                       << "Ref.   = " << test.ref_cs << "\n\t"
                       << "CepGen = " << new_cs << "\n\t"
                       << "Ratio: " << ratio << "\n\t"
                       << "Pull: " << pull << ".\n\t"
                       << "Computation time: " << tmr.elapsed() * 1.e3 << " ms.";
      tmr.reset();

      const string test_res = cepgen::utils::format(
          "%-40s\tref=%g\tgot=%g\tratio=%g\tpull=%+10.5f", test.filename.c_str(), test.ref_cs, new_cs, ratio, pull);
      CG_TEST(fabs(pull) < num_sigma, filename);
      gen.runParameters().clearProcess();
    } catch (const cepgen::utils::RunAbortedException&) {
      CG_TEST_SUMMARY;
    } catch (const cepgen::Exception& e) {
      CG_LOG << "Test \"" << test.filename << "\" (located at " << filename << ") failed.";
      cepgen::Exception(e).dump();
    }
  }
  CG_TEST_SUMMARY;
}
