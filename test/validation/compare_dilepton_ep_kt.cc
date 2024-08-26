/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"
#include "Comparator.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  int num_gen;
  vector<string> processes;
  string filename, plotter;
  bool ratio_plot;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument(
          "processes,P", "processes to generate", &processes, vector<string>{"lpair", "pptoff", "mg5_aMC"})
      .addOptionalArgument("num-gen,n", "number of events to generate", &num_gen, 10'000)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "root")
      .addOptionalArgument("ratio,r", "draw the ratio plot", &ratio_plot, false)
      .addOptionalArgument("filename,f", "output base filename", &filename, "validation/comparison_dilepton_ep_kt_")
      .parse();

  struct comparison_t : public cepgen::validation::Comparator {
    using cepgen::validation::Comparator::Comparator;
    void initialise() override {
      (*this)
          .book("invmass", "$m(l^{+}l^{-})$", "GeV", cepgen::utils::Hist1D(50, {10., 160.}))
          .book("ptpair", "$p_{T}(l^{+}l^{-})$", "GeV", cepgen::utils::Hist1D(50, {0., 5.}))
          .book("ptlead", "$p_{T}^{lead}$", "GeV", cepgen::utils::Hist1D(50, {0., 50.}))
          .book("ptsublead", "$p_{T}^{sublead}$", "GeV", cepgen::utils::Hist1D(50, {0., 50.}))
          .book("etalead", "$\\eta^{lead}$", "", cepgen::utils::Hist1D(50, {-2.5, 2.5}))
          .book("etasublead", "$\\eta^{sublead}$", "", cepgen::utils::Hist1D(50, {-2.5, 2.5}))
          .book("acop", "$1-|\\Delta\\phi(l^{+}l^{-})/\\pi|$", "", cepgen::utils::Hist1D(50, {0., 0.5}))
          .book("mx", "$M_{X}$", "GeV", cepgen::utils::Hist1D(50, {0., 1000.}));
      for (const auto& plot : {"invmass", "ptpair", "ptlead", "ptsublead", "acop"})
        drawMode(plot) |= cepgen::utils::Drawer::Mode::logy;
    }
    void process(const cepgen::Event& evt) override {
      const auto &cm = evt(cepgen::Particle::Role::Intermediate).at(0).momentum(),
                 &pl1 = evt(cepgen::Particle::Role::CentralSystem).at(0).momentum(),
                 &pl2 = evt(cepgen::Particle::Role::CentralSystem).at(1).momentum();
      cepgen::Momentum pl_lead, pl_sublead;
      if (pl1.pt() > pl2.pt())
        pl_lead = pl1, pl_sublead = pl2;
      else
        pl_lead = pl2, pl_sublead = pl1;
      (*this)
          .fill("invmass", cm.mass())
          .fill("ptpair", cm.pt())
          .fill("ptlead", pl_lead.pt())
          .fill("etalead", pl_lead.eta())
          .fill("ptsublead", pl_sublead.pt())
          .fill("etasublead", pl_sublead.eta())
          .fill("acop", 1. - fabs(pl1.deltaPhi(pl2) * M_1_PI))
          .fill("mx", evt(cepgen::Particle::Role::OutgoingBeam1).at(0).momentum().mass());
    }
  };

  cepgen::Generator gen;
  auto comp = comparison_t(
      gen,
      cepgen::ParametersList()
          .set("topLabel", "SD $\\gamma\\gamma \\rightarrow l^{+}l^{-}$ (13.6 TeV), $p_{T}^{l} > 10$ GeV, $k_{T}$"s)
          .set("numEvents", num_gen)
          .set("pathTemplate", filename)
          .set("plotter", cepgen::ParametersList().setName(plotter).feed(plotter).set("format", "png,pdf"s)));

  auto& pars = gen.runParameters();
  for (const auto& proc_name : processes) {
    auto proc = proc_name;
    if (!cepgen::utils::contains(cepgen::ProcessFactory::get().modules(), proc_name))
      continue;
    if (proc_name == "mg5_aMC")
      proc += "<process:'a a > mu- mu+'";
    pars.setProcess(cepgen::ProcessFactory::get().build(
        proc,
        cepgen::ParametersList().set(
            "kinematicsGenerator",
            cepgen::PhaseSpaceGeneratorFactory::get().describeParameters("kt:2to4").parameters())));
    pars.process().kinematics().setParameters(cepgen::ParametersList()
                                                  .set<vector<int> >("pdgIds", {2212, 11})
                                                  .set<vector<double> >("pz", {7000., 50.})
                                                  .set<int>("mode", 1 /* elastic-elastic */)
                                                  .set<double>("ptmin", 2.5));
    comp.loop(proc_name);
  }
  return 0;
}
