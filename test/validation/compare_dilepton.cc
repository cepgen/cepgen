/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2026  Laurent Forthomme
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
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/Message.h"
#include "Comparator.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  int num_gen, num_threads, mode;
  double sqrt_s, pt_min;
  vector<string> processes;
  string filename, plotter;
  bool coll_factor, ratio_plot;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("coll,c", "use collinear factorisation?", &coll_factor, false)
      .addOptionalArgument("processes,P", "processes to generate", &processes, vector{"lpair"s, "pptoff"s, "mg5_aMC"s})
      .addOptionalArgument("num-gen,n", "number of events to generate", &num_gen, 50'000)
      .addOptionalArgument("num-threads,t", "number of threads to generate events", &num_threads, 4)
      .addOptionalArgument("mode", "kinematics mode (1=el-el, 2=el-inel, 3=inel-el, 4=inel-inel)", &mode, 3)
      .addOptionalArgument("sqrts,s", "centre of mass energy (in GeV)", &sqrt_s, 14.e3)
      .addOptionalArgument("ptmin", "minimum slepton pT (in GeV)", &pt_min, 10.)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "root")
      .addOptionalArgument("ratio,r", "draw the ratio plot", &ratio_plot, false)
      .addOptionalArgument("filename,f", "output base filename", &filename, "validation/comparison_dilepton_")
      .parse();

  std::string ps_gen_name, ps_gen_label;
  if (coll_factor) {
    ps_gen_name = "coll";
    ps_gen_label = "coll.";
  } else {
    ps_gen_name = "kt";
    ps_gen_label = "$k_{T}$";
  }

  struct comparison_t final : cepgen::validation::Comparator {
    explicit comparison_t(cepgen::Generator& gen, const cepgen::ParametersList& params, bool _coll_factor)
        : Comparator(gen, params), coll_factor(_coll_factor) {}
    void initialise() override {
      (*this)
          .book("invmass", "$m(l^{+}l^{-})$", "GeV", cepgen::utils::Hist1D(50, {10., 510.}))
          .book("ptpair", "$p_{T}(l^{+}l^{-})$", "GeV", cepgen::utils::Hist1D(50, {0., coll_factor ? 5. : 50.}))
          .book("ptlead", "$p_{T}^{lead}$", "GeV", cepgen::utils::Hist1D(50, {0., 100.}))
          .book("ptsublead", "$p_{T}^{sublead}$", "GeV", cepgen::utils::Hist1D(50, {0., 100.}))
          .book("etalead", "$\\eta^{lead}$", "", cepgen::utils::Hist1D(50, {-2.5, 2.5}))
          .book("etasublead", "$\\eta^{sublead}$", "", cepgen::utils::Hist1D(50, {-2.5, 2.5}))
          .book("acop", "$1-|\\Delta\\phi(l^{+}l^{-})/\\pi|$", "", cepgen::utils::Hist1D(50, {0., 1.}))
          .book("mx", "$M_{X}$", "GeV", cepgen::utils::Hist1D(50, {0., 1000.}));
      for (const auto& plot : {"invmass", "ptpair", "ptlead", "ptsublead", "acop", "mx"})
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
    bool coll_factor;
  };

  cepgen::Generator gen;
  auto comp = comparison_t(
      gen,
      cepgen::ParametersList()
          .set("topLabel",
               cepgen::validation::Comparator::modeFromId(mode) +
                   cepgen::utils::format(" $\\gamma\\gamma \\rightarrow l^{+}l^{-}$ (%g TeV)", sqrt_s * 1.e-3) +
                   (pt_min > 0. ? cepgen::utils::format(", $p_{T}^{l} >$ %g GeV, "s, pt_min) : ""s) + ps_gen_label)
          .set("numEvents", num_gen)
          .set("numThreads", num_threads)
          .set("pathTemplate", filename)
          .set("plotter", cepgen::ParametersList().setName(plotter).feed(plotter).set("format", "png,pdf"s)),
      coll_factor);

  auto& pars = gen.runParameters();
  for (const auto& proc_name : processes) {
    const auto module_name = cepgen::utils::split(proc_name, '<').at(0);
    auto proc = proc_name;
    if (!cepgen::utils::contains(cepgen::ProcessFactory::get().modules(), module_name)) {
      CG_WARNING("compare_dilepton") << "Process module '" << module_name << "' is not registered; skipping.";
      continue;
    }
    if (module_name == "mg5_aMC")
      proc += "<process:'a a > mu- mu+'";
    pars.setProcess(
        cepgen::ProcessFactory::get().build(proc, cepgen::ParametersList().set("partonsGenerator", ps_gen_name)));
    pars.process().kinematics().setParameters(cepgen::ParametersList()
                                                  .set<vector<int> >("pdgIds", {2212, 2212})
                                                  .set<double>("sqrtS", sqrt_s)
                                                  .set<int>("mode", mode)
                                                  .set<double>("ptmin", pt_min));
    comp.loop(proc_name);
  }
  return 0;
}
