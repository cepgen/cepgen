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

#include <TGraph.h>
#include <TMultiGraph.h>

#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

using namespace std;

int main(int argc, char* argv[]) {
  double qmin, qmax;
  int num_points;
  string output_file;
  bool logy, draw_grid;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("qmin,m", "minimum virtuality (GeV)", &qmin, 1.)
      .addOptionalArgument("qmax,M", "maximum virtuality (GeV)", &qmax, 101.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 100)
      .addOptionalArgument("output,o", "output file name", &output_file, "alphas.scan.output.txt")
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .parse();

  cepgen::initialise();

  struct alpha_t {
    string name;
    vector<double> vals;
    TGraph graph;
  };
  vector<alpha_t> alphas, alphaem;

  vector<double> qvals(num_points);
  for (int i = 0; i < num_points; ++i)
    qvals[i] = qmin + (qmax - qmin) * i / num_points;

  // alphaS(Q) modellings part
  size_t i = 0;
  for (const auto& mod : cepgen::AlphaSFactory::get().modules()) {
    const auto& algo = cepgen::AlphaSFactory::get().build(
        mod /*, cepgen::ParametersList().set<double>("asmur", 0.35).set<double>("mur", 1.4142)*/);
    TGraph graph;
    graph.SetName(mod.c_str());
    alphas.emplace_back(alpha_t{mod, vector<double>(num_points), graph});
    auto& as = alphas[i++];
    for (size_t j = 0; j < qvals.size(); ++j) {
      const auto val = (*algo)(qvals[j]);
      as.vals[j] = val;
      as.graph.SetPoint(j, qvals[j], val);
    }
  }
  // alphaEM(Q) modellings part
  i = 0;
  for (const auto& mod : cepgen::AlphaEMFactory::get().modules()) {
    const auto& algo = cepgen::AlphaEMFactory::get().build(mod);
    TGraph graph;
    graph.SetName(mod.c_str());
    alphaem.emplace_back(alpha_t{mod, vector<double>(num_points), graph});
    auto& aem = alphaem[i++];
    for (size_t j = 0; j < qvals.size(); ++j) {
      const auto val = (*algo)(qvals[j]);
      aem.vals[j] = val;
      aem.graph.SetPoint(j, qvals[j], val);
    }
  }

  // output ascii file
  ofstream out(output_file);
  out << "#";
  for (const auto& smp : alphas)
    out << "\t" << smp.name;
  for (const auto& smp : alphaem)
    out << "\t" << smp.name;
  for (size_t i = 0; i < qvals.size(); ++i) {
    out << "\n" << qvals[i];
    for (const auto& smp : alphas)
      out << "\t" << smp.vals[i];
    for (const auto& smp : alphaem)
      out << "\t" << smp.vals[i];
  }

  // drawing part
  const auto top_label = cepgen::utils::s("CepGen #alpha_{S,EM} modelling", alphas.size() + alphaem.size(), false);
  cepgen::ROOTCanvas c("comp_alphas_alphaem", top_label.c_str());
  c.SetLegendX1(0.15);
  if (draw_grid)
    c.SetGrid(true, true);
  TMultiGraph mg;
  vector<TH1*> numers(alphas.size() + alphaem.size());
  for (size_t i = 0; i < alphas.size(); ++i) {
    auto& graph = alphas[i].graph;
    graph.SetLineColor(cepgen::ROOTCanvas::colours[i]);
    mg.Add(&graph);
    numers[i] = graph.GetHistogram();
    const auto descr = cepgen::utils::replace_all(cepgen::AlphaSFactory::get().describe(alphas[i].name),
                                                  {{" alphaS", ""}, {" evolution algorithm", ""}});
    c.AddLegendEntry(&graph, descr.c_str(), "l");
  }
  for (size_t i = 0; i < alphaem.size(); ++i) {
    auto& graph = alphaem[i].graph;
    graph.SetLineColor(cepgen::ROOTCanvas::colours[i]);
    graph.SetLineStyle(2);
    mg.Add(&graph);
    numers[i] = graph.GetHistogram();
    const auto descr = cepgen::utils::replace_all(cepgen::AlphaEMFactory::get().describe(alphaem[i].name),
                                                  {{" alphaS", ""}, {" evolution algorithm", ""}});
    c.AddLegendEntry(&graph, descr.c_str(), "l");
  }
  mg.Draw("al");
  mg.GetHistogram()->SetTitle(";Q (GeV);#alpha_{S,EM}(Q)");
  mg.GetXaxis()->SetRangeUser(*qvals.begin(), *qvals.rbegin());
  c.Prettify(mg.GetHistogram());
  c.SetLogx();
  if (logy) {
    c.SetLogy();
    mg.SetMinimum(1.e-3);
  }
  c.Save("pdf");

  return 0;
}
