#include "CepGen/Generator.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Physics/AlphaS.h"
#include "CepGen/Utils/ArgumentsParser.h"

#include <fstream>

#include "Canvas.h"
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

  vector<double> qvals(num_points);
  vector<pair<string, vector<double> > > alphas_vals;
  for (int i = 0; i < num_points; ++i)
    qvals[i] = qmin + (qmax - qmin) * i / num_points;

  size_t i = 0;
  vector<TGraph> v_graphs;
  for (const auto& mod : cepgen::AlphaSFactory::get().modules()) {
    const auto& algo = cepgen::AlphaSFactory::get().build(mod);
    alphas_vals.emplace_back(make_pair(mod, vector<double>(num_points)));
    v_graphs.emplace_back();
    for (size_t j = 0; j < qvals.size(); ++j) {
      const auto val = (*algo)(qvals[j]);
      alphas_vals[i].second[j] = val;
      v_graphs[i].SetPoint(j, qvals[j], val);
    }
    ++i;
  }

  ofstream out(output_file);
  out << "#";
  for (const auto& smp : alphas_vals)
    out << "\t" << smp.first;
  for (size_t i = 0; i < qvals.size(); ++i) {
    out << "\n" << qvals[i];
    for (const auto& smp : alphas_vals)
      out << "\t" << smp.second[i];
  }

  cepgen::Canvas c("test_alphas");
  c.SetLegendX1(0.15);
  TMultiGraph mg;
  for (size_t i = 0; i < alphas_vals.size(); ++i) {
    v_graphs[i].SetLineColor(cepgen::Canvas::colours[i]);
    mg.Add(&v_graphs[i]);
    c.AddLegendEntry(&v_graphs[i], cepgen::AlphaSFactory::get().describe(alphas_vals[i].first).c_str());
  }
  mg.Draw("al");
  mg.GetHistogram()->SetTitle(";Q (GeV);#alpha_{S}(Q)");
  c.Prettify(mg.GetHistogram());
  c.SetLogx();
  c.Save("pdf");

  return 0;
}
