/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#include "CepGen/Core/EventExporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/Message.h"
#include "CepGenAddOns/ROOTWrapper/ROOTTreeInfo.h"

using namespace std;

int main(int argc, char* argv[]) {
  bool keep_file;
  string proc_name, tmp_filename;
  int num_gen;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("keep-file,k", "keep the output TTree", &keep_file, false)
      .addOptionalArgument("process,p", "process to generate", &proc_name, "lpair")
      .addOptionalArgument("filename,f", "temporary filename", &tmp_filename, "/tmp/cepgen_test.root")
      .addOptionalArgument("num-gen,n", "number of events to generate", &num_gen, 10)
      .parse();

  double cross_sec, cross_sec_unc;
  {  // generation + tree building part
    cepgen::Generator gen;
    auto& pars = gen.parametersRef();
    pars.setProcess(cepgen::proc::ProcessFactory::get().build(proc_name));
    pars.process().setKinematics(cepgen::Kinematics(cepgen::ParametersList()
                                                        .set<vector<int> >("pdgIds", {2212, 2212})
                                                        .set<double>("sqrtS", 13.e3)
                                                        .set<int>("mode", 1)
                                                        .set<double>("ptmin", 25.)));
    pars.addEventExporter(cepgen::EventExporterFactory::get().build(
        "root_tree", cepgen::ParametersList().set<string>("filename", tmp_filename)));
    CG_LOG << &pars;
    gen.generate(num_gen);
    cross_sec = gen.crossSection();
    cross_sec_unc = gen.crossSectionError();
  }

  {  // tree analysis part
    auto* file = TFile::Open(tmp_filename.c_str());
    ROOT::CepGenRun run_info;
    run_info.attach(file);
    if (run_info.xsect != cross_sec) {
      CG_LOG << "Invalid cross section retrieved from run tree: " << run_info.xsect << " != " << cross_sec << ".";
      return -1;
    }
    if (run_info.errxsect != cross_sec_unc) {
      CG_LOG << "Invalid cross section uncertainty retrieved from run tree: " << run_info.errxsect
             << " != " << cross_sec_unc << ".";
      return -1;
    }
    ROOT::CepGenEvent evt_info;
    evt_info.attach(file);
    if (!evt_info.tree()) {
      CG_LOG << "Failed to retrieve an events tree.";
      return -1;
    }
    const auto num_gen_eff = evt_info.tree()->GetEntriesFast();
    if (num_gen_eff != num_gen) {
      CG_LOG << "Invalid number of events generated: " << num_gen_eff << " != " << num_gen << ".";
      return -1;
    }
  }

  {  // tree removal part
    if (!keep_file && !fs::remove(tmp_filename)) {
      CG_LOG << "Failed to remove the temporary file \"" << tmp_filename << "\".";
      return -1;
    }
  }

  return 0;
}
