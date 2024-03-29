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
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  int num_gen;
  vector<string> processes;
  string filename, plotter;
  bool ratio_plot;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("processes,p", "processes to generate", &processes, vector<string>{"pptoff", "mg5_aMC"})
      .addOptionalArgument("num-gen,n", "number of events to generate", &num_gen, 10'000)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "root")
      .addOptionalArgument("ratio,r", "draw the ratio plot", &ratio_plot, false)
      .addOptionalArgument(
          "filename,f",
          "output base filename",
          &filename,
          fs::path(cepgen::utils::env::get("CEPGEN_PATH", ".")) / "validation" / "comparison_dilepton_coll_")
      .parse();

  vector<cepgen::utils::Hist1D> h_invmass(processes.size(),
                                          cepgen::utils::Hist1D(50, {10., 510.}, "", "$m(l^{+}l^{-})$ (GeV)")),
      h_ptpair(processes.size(), cepgen::utils::Hist1D(50, {0., 5.}, "", "$p_{T}(l^{+}l^{-})$ (GeV)")),
      h_ptlead(processes.size(), cepgen::utils::Hist1D(50, {0., 100.}, "", "$p_{T}^{lead}$ (GeV)")),
      h_ptsublead(processes.size(), cepgen::utils::Hist1D(50, {0., 100.}, "", "$p_{T}^{sublead}$ (GeV)")),
      h_etalead(processes.size(), cepgen::utils::Hist1D(50, {-2.5, 2.5}, "", "$\\eta^{lead}$")),
      h_etasublead(processes.size(), cepgen::utils::Hist1D(50, {-2.5, 2.5}, "", "$\\eta^{sublead}$")),
      h_acop(processes.size(), cepgen::utils::Hist1D(50, {0., 1.}, "", "1-|\\Delta\\phi(l^{+}l^{-})/\\pi|")),
      h_mx(processes.size(), cepgen::utils::Hist1D(50, {0., 1000.}, "", "M_{X} (GeV)"));

  cepgen::Generator gen;
  auto& pars = gen.runParameters();
  size_t i = 0;
  const string plot_title = "SD $\\gamma\\gamma \\rightarrow l^{+}l^{-}$ (13.6 TeV), $p_{T}^{l} > 10$ GeV, coll.";
  for (const auto& proc_name : processes) {
    auto proc = proc_name;
    if (proc_name == "mg5_aMC")
      proc += "<process:'a a > mu- mu+'";
    else if (proc_name == "pptoff")
      proc += "<method:0";
    pars.setProcess(cepgen::ProcessFactory::get().build(
        proc,
        cepgen::ParametersList().set(
            "kinematicsGenerator",
            cepgen::PhaseSpaceGeneratorFactory::get().describeParameters("coll2to4").parameters())));
    pars.process().kinematics().setParameters(cepgen::ParametersList()
                                                  .set<vector<int> >("pdgIds", {2212, 2212})
                                                  .set<double>("sqrtS", 13.6e3)
                                                  .set<int>("mode", 3 /* inelastic-elastic */)
                                                  .set<double>("ptmin", 10.));
    const auto cs = gen.computeXsection();
    CG_LOG << "Cross section computed for process '" << proc_name << "': " << cs << " pb.";
    const auto weight = (double)cs / num_gen;
    gen.generate(num_gen, [&](const cepgen::Event& evt, size_t) {
      const auto &cm = evt(cepgen::Particle::Role::Intermediate).at(0).momentum(),
                 &px = evt(cepgen::Particle::Role::OutgoingBeam1).at(0).momentum(),
                 &pl1 = evt(cepgen::Particle::Role::CentralSystem).at(0).momentum(),
                 &pl2 = evt(cepgen::Particle::Role::CentralSystem).at(1).momentum();
      h_invmass[i].fill(cm.mass(), weight);
      h_ptpair[i].fill(cm.pt(), weight);
      cepgen::Momentum pl_lead, pl_sublead;
      if (pl1.pt() > pl2.pt())
        pl_lead = pl1, pl_sublead = pl2;
      else
        pl_lead = pl2, pl_sublead = pl1;
      h_ptlead[i].fill(pl_lead.pt(), weight);
      h_etalead[i].fill(pl_lead.eta(), weight);
      h_ptsublead[i].fill(pl_sublead.pt(), weight);
      h_etasublead[i].fill(pl_sublead.eta(), weight);
      h_acop[i].fill(1. - fabs(pl1.deltaPhi(pl2) * M_1_PI), weight);
      h_mx[i].fill(px.mass());
    });
    ++i;
  }
  if (!plotter.empty()) {
    auto plt = cepgen::DrawerFactory::get().build(plotter, cepgen::ParametersList().set<string>("format", "png,pdf"));
    cepgen::utils::Drawer::Mode dm = cepgen::utils::Drawer::Mode::nostack | cepgen::utils::Drawer::Mode::grid;
    if (ratio_plot)
      dm |= cepgen::utils::Drawer::Mode::ratio;

    for (auto& plot : vector<pair<string, vector<cepgen::utils::Hist1D> > >{{"invmass", h_invmass},
                                                                            {"ptpair", h_invmass},
                                                                            {"ptlead", h_ptlead},
                                                                            {"etalead", h_etalead},
                                                                            {"ptsublead", h_ptsublead},
                                                                            {"etasublead", h_etasublead},
                                                                            {"acop", h_acop},
                                                                            {"mx", h_mx}}) {
      cepgen::utils::DrawableColl coll;
      auto mode = dm;
      if (plot.first == "etalead" || plot.first == "etasublead") {
      } else
        mode |= cepgen::utils::Drawer::Mode::logy;
      size_t i = 0;
      for (auto& gr : plot.second) {
        gr.xAxis().setLabel(gr.title());
        gr.yAxis().setLabel("d$\\sigma/dx");
        string chi2_info;
        if (i > 0) {  // take first plot as reference for chi^2 test
          size_t ndf;
          const auto chi2 = gr.chi2test(plot.second.at(0), ndf);
          chi2_info = cepgen::utils::format(", $\\chi^{2}$/ndf = %.2g/%zu", chi2, ndf);
        }
        gr.setTitle(processes.at(i) + chi2_info);
        coll.emplace_back(&gr);
        ++i;
      }
      plt->draw(coll, filename + plot.first, plot_title, mode);
    }
  }
  return 0;
}
