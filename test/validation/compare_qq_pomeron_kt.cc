/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "Comparator.h"

using namespace std;
using namespace std::string_literals;

/// Compare the elastic gamma-gamma, gamma-pomeron and pomeron-pomeron production of a quark pair over the same phase space
int main(int argc, char* argv[]) {
  int num_gen, pair;
  string filename, plotter;
  double sqrts, ptmin, qtmax;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("pair,q", "PDG id of the produced quark flavour", &pair, 4 /* charm */)
      .addOptionalArgument("num-gen,n", "number of events to generate", &num_gen, 50'000)
      .addOptionalArgument("sqrts,s", "centre of mass energy (GeV)", &sqrts, 13.6e3)
      .addOptionalArgument("ptmin", "minimal quark transverse momentum (GeV)", &ptmin, 3.)
      .addOptionalArgument("qtmax", "maximal parton transverse virtuality (GeV)", &qtmax, 50.)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "root")
      .addOptionalArgument("filename,f", "output base filename", &filename, "validation/comparison_qq_pomeron_kt_")
      .parse();

  struct comparison_t : public cepgen::validation::Comparator {
    using cepgen::validation::Comparator::Comparator;
    void initialise() override {
      (*this)
          .book("invmass", "$m(q\\bar{q})$", "GeV", cepgen::utils::Hist1D(50, {0., 100.}))
          .book("ptpair", "$p_{T}(q\\bar{q})$", "GeV", cepgen::utils::Hist1D(50, {0., 50.}))
          .book("ptlead", "$p_{T}^{lead}$", "GeV", cepgen::utils::Hist1D(50, {0., 50.}))
          .book("ptsublead", "$p_{T}^{sublead}$", "GeV", cepgen::utils::Hist1D(50, {0., 50.}))
          .book("etalead", "$\\eta^{lead}$", "", cepgen::utils::Hist1D(50, {-2.5, 2.5}))
          .book("etasublead", "$\\eta^{sublead}$", "", cepgen::utils::Hist1D(50, {-2.5, 2.5}))
          .book("acop", "$1-|\\Delta\\phi(q\\bar{q})/\\pi|$", "", cepgen::utils::Hist1D(50, {0., 0.5}));
      for (const auto& plot : {"invmass", "ptpair", "ptlead", "ptsublead", "acop"})
        drawMode(plot) |= cepgen::utils::Drawer::Mode::logy;
    }
    void process(const cepgen::Event& evt) override {
      const auto &cm = evt(cepgen::Particle::Role::Intermediate).at(0).momentum(),
                 &pq1 = evt(cepgen::Particle::Role::CentralSystem).at(0).momentum(),
                 &pq2 = evt(cepgen::Particle::Role::CentralSystem).at(1).momentum();
      const auto &pq_lead = pq1.pt() > pq2.pt() ? pq1 : pq2, &pq_sublead = pq1.pt() > pq2.pt() ? pq2 : pq1;
      (*this)
          .fill("invmass", cm.mass())
          .fill("ptpair", cm.pt())
          .fill("ptlead", pq_lead.pt())
          .fill("etalead", pq_lead.eta())
          .fill("ptsublead", pq_sublead.pt())
          .fill("etasublead", pq_sublead.eta())
          .fill("acop", 1. - fabs(pq1.deltaPhi(pq2) * M_1_PI));
    }
  };

  cepgen::Generator gen;

  std::ostringstream os_pdg;
  os_pdg << cepgen::PDG::get().name(pair);
  auto comp = comparison_t(
      gen,
      cepgen::ParametersList()
          .set("topLabel",
               cepgen::utils::format("El. $%s\\bar{%s}$ prod. (%g TeV), $p_{T}^{q} > %g$ GeV, $k_{T}$",
                                     os_pdg.str().data(),
                                     os_pdg.str().data(),
                                     sqrts * 1.e-3,
                                     ptmin))
          .set("numEvents", num_gen)
          .set("pathTemplate", filename)
          .set("plotter", cepgen::ParametersList().setName(plotter).feed(plotter).set("format", "png,pdf"s)));

  // incoming parton types, defined for each beam
  const vector<std::pair<string, vector<string> > > samples{{"$\\gamma\\gamma$", {"BudnevElastic", "BudnevElastic"}},
                                                            {"$\\gamma IP$", {"BudnevElastic", "PomeronElastic"}},
                                                            {"$IPIP$", {"PomeronElastic", "PomeronElastic"}}};

  auto& pars = gen.runParameters();
  for (const auto& [sample_name, fluxes] : samples) {
    pars.setProcess(cepgen::ProcessFactory::get().build(
        "pptoff", cepgen::ParametersList().set<string>("partonsGenerator", "kt").set<int>("pair", pair)));
    pars.process().kinematics().setParameters(cepgen::ParametersList()
                                                  .set<vector<int> >("pdgIds", {2212, 2212})
                                                  .set<double>("sqrtS", sqrts)
                                                  .set<int>("mode", 1 /* elastic-elastic */)
                                                  .set<vector<string> >("partonFluxes", fluxes)
                                                  .set<double>("ptmin", ptmin)
                                                  .set<cepgen::Limits>("qt", cepgen::Limits{0., qtmax}));
    comp.loop(sample_name);
  }
  return 0;
}
