/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2026  Laurent Forthomme
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

#include <CepGen/Core/RunParameters.h>
#include <CepGen/EventFilter/EventExporter.h>
#include <CepGen/Generator.h>
#include <CepGen/Modules/EventExporterFactory.h>
#include <CepGen/Modules/ProcessFactory.h>
#include <CepGen/Process/Process.h>
#include <CepGen/Utils/ArgumentsParser.h>
#include <CepGen/Utils/Filesystem.h>
#include <CepGen/Utils/Message.h>
#include <CepGen/Utils/Test.h>

#include "CepGenRoot/ROOTTreeInfo.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  bool keep_file;
  string proc_name, tmp_filename;
  int num_gen;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("keep-file,k", "keep the output TTree", &keep_file, false)
      .addOptionalArgument("process,p", "process to generate", &proc_name, "lpair")
      .addOptionalArgument("filename,f", "temporary filename", &tmp_filename, "cepgen_test.root")
      .addOptionalArgument("num-gen,n", "number of events to generate", &num_gen, 10)
      .parse();

  cepgen::utils::TemporaryDirectory temp_dir("cepgen_root_test", keep_file);
  const auto full_path = temp_dir.path() / tmp_filename;
  double cross_sec, cross_sec_unc;
  {  // generation + tree building part
    cepgen::Generator gen;
    auto& pars = gen.runParameters();
    pars.setProcess(cepgen::ProcessFactory::get().build(proc_name));
    pars.process().kinematics().setParameters(cepgen::ParametersList()
                                                  .set("pdgIds"s, std::vector{2212, 2212})
                                                  .set("sqrtS"s, 13.6e3)
                                                  .set("mode"s, 1)
                                                  .set("ptmin"s, 25.));
    pars.addEventExporter(cepgen::EventExporterFactory::get().build(
        "root_tree"s, cepgen::ParametersList().set("filename", full_path.string())));
    gen.generate(num_gen);
    cross_sec = gen.crossSection();
    cross_sec_unc = gen.crossSectionError();
  }

  {  // tree analysis part
    TFile file(full_path.string().c_str());
    ROOT::CepGenRun run_info;
    run_info.attach(&file);

    CG_TEST_EQUIV(run_info.xsect, cross_sec, "cross section from run tree");
    CG_TEST_EQUIV(run_info.errxsect, cross_sec_unc, "cross section uncertainty from run tree");

    ROOT::CepGenEvent evt_info;
    evt_info.attach(&file);
    CG_TEST(evt_info.tree(), "events tree present");
    if (evt_info.tree())
      CG_TEST(evt_info.tree()->GetEntriesFast() == num_gen, "number of events generated");
  }

  CG_TEST_SUMMARY;
}
