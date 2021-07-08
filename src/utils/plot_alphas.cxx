#include "CepGen/Generator.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Physics/AlphaS.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/ROOTWrapper/Canvas.h"

#include <fstream>

#include <TMultiGraph.h>
#include <TGraph.h>

using namespace std;

int main(int argc, char* argv[]) {
  double qmin, qmax;
  int num_points;
  string output_file;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("qmin,m", "minimum virtuality (GeV)", &qmin, 1.)
      .addOptionalArgument("qmax,M", "maximum virtuality (GeV)", &qmax, 101.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 100)
      .addOptionalArgument("output,o", "output file name", &output_file, "alphas.scan.output.txt")
      .parse();

  cepgen::initialise();

  struct alphas_t {
    string name;
    vector<double> vals;
    TGraph graph;
  };
  vector<alphas_t> alphas;

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
    alphas.emplace_back(alphas_t{mod, vector<double>(num_points), graph});
    auto& as = alphas[i++];
    for (size_t j = 0; j < qvals.size(); ++j) {
      const auto val = (*algo)(qvals[j]);
      as.vals[j] = val;
      as.graph.SetPoint(j, qvals[j], val);
    }
  }

  // output ascii file
  ofstream out(output_file);
  out << "#";
  for (const auto& smp : alphas)
    out << "\t" << smp.name;
  for (size_t i = 0; i < qvals.size(); ++i) {
    out << "\n" << qvals[i];
    for (const auto& smp : alphas)
      out << "\t" << smp.vals[i];
  }

  // drawing part
  const auto top_label = cepgen::utils::s("CepGen #alpha_{S} modelling", alphas.size(), false);
  cepgen::Canvas c("test_alphas", top_label.c_str(), alphas.size() > 1);
  c.SetLegendX1(0.15);
  TMultiGraph mg;
  vector<TH1*> numers(alphas.size());
  for (size_t i = 0; i < alphas.size(); ++i) {
    auto& graph = alphas[i].graph;
    graph.SetLineColor(cepgen::Canvas::colours[i]);
    mg.Add(&graph);
    numers[i] = graph.GetHistogram();
    auto descr = cepgen::AlphaSFactory::get().describe(alphas[i].name);
    cepgen::utils::replace_all(descr, " alphaS", "");
    cepgen::utils::replace_all(descr, " evolution algorithm", "");
    c.AddLegendEntry(&graph, descr.c_str(), "l");
  }
  c.RatioPlot(alphas[0].graph.GetHistogram(), numers);
  mg.Draw("al");
  mg.GetHistogram()->SetTitle(";Q (GeV);#alpha_{S}(Q)");
  mg.GetXaxis()->SetRangeUser(*qvals.begin(), *qvals.rbegin());
  c.Prettify(mg.GetHistogram());
  c.SetLogx();
  c.Save("pdf");

  return 0;
}
