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
#include "CepGen/Utils/Test.h"
#include "CepGen/Utils/Timer.h"

using namespace std;

int main(int argc, char* argv[]) {
  double num_sigma;
  string cfg_filename;
  string integrator;
  bool verbose;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("cfg,f", "configuration file", &cfg_filename, "test/physics/test_processes.cfg")
      .addOptionalArgument("verbose,v", "verbose mode", &verbose, false)
      .addOptionalArgument("num-sigma,n", "max. number of std.dev.", &num_sigma, 3.)
      .addOptionalArgument("integrator,i", "type of integrator used", &integrator, "Vegas")
      .parse();

  if (!verbose)
    CG_LOG_LEVEL(warning);

  cepgen::utils::Timer tmr;
  cepgen::Generator gen;

  CG_LOG << "Testing with " << integrator << " integrator.";

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

  for (const auto& test : tests) {
    const std::string filename = "test/physics/test_processes/" + test.filename + "_cfg.py";
    try {
      gen.parametersRef().clearProcess();
      gen.setParameters(cepgen::card::Handler::parse(filename));

      gen.parametersRef().par_integrator.setName<std::string>(integrator);
      CG_DEBUG("main") << "Process: " << gen.parameters()->processName() << "\n\t"
                       << "File: " << filename << "\n\t"
                       << "Configuration time: " << tmr.elapsed() * 1.e3 << " ms.";

      tmr.reset();

      double new_cs, err_new_cs;
      gen.computeXsection(new_cs, err_new_cs);

      const double ratio = new_cs / test.ref_cs;
      const double err_ratio = ratio * hypot(err_new_cs / new_cs, test.err_ref_cs / test.ref_cs);
      const double pull = (new_cs - test.ref_cs) / hypot(err_new_cs, test.err_ref_cs);

      CG_DEBUG("main") << "Computed cross section:\n\t"
                       << "Ref.   = " << test.ref_cs << " +/- " << test.err_ref_cs << "\n\t"
                       << "CepGen = " << new_cs << " +/- " << err_new_cs << "\n\t"
                       << "Ratio: " << ratio << " +/- " << err_ratio << "\n\t"
                       << "Pull: " << pull << ".\n\t"
                       << "Computation time: " << tmr.elapsed() * 1.e3 << " ms.";
      tmr.reset();

      const string test_res = cepgen::utils::format(
          "%-40s\tref=%g\tgot=%g\tratio=%g\tpull=%+10.5f", test.filename.c_str(), test.ref_cs, new_cs, ratio, pull);
      CG_TEST(fabs(pull) < num_sigma, filename);
    } catch (const cepgen::utils::RunAbortedException&) {
      CG_TEST_SUMMARY;
    } catch (const cepgen::Exception& e) {
      CG_LOG << "Test \"" << test.filename << "\" (located at " << filename << ") failed.";
      cepgen::Exception(e).dump();
    }
  }
  CG_TEST_SUMMARY;
}
