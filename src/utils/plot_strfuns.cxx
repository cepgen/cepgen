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
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<int> strfun_types;
  double q2, xmin, xmax;
  int num_points;
  string output_file;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("sf,s", "structure functions modelling", &strfun_types)
      .addOptionalArgument("q2,q", "parton virtuality (GeV^2)", &q2, 10.)
      .addOptionalArgument("xmax,m", "minimal Bjorken x", &xmin, 1.e-7)
      .addOptionalArgument("xmax,M", "maximal Bjorken x", &xmax, 1.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "strfuns.scan.output.txt")
      .parse();

  cepgen::initialise();

  ofstream out(output_file);
  out << "# structure functions: ";
  string sep;
  for (const auto& sf_type : strfun_types)
    out << sep << sf_type, sep = ", ";
  out << "\n"
      << "# x in [" << xmin << ", " << xmax << "]\n";

  vector<unique_ptr<cepgen::strfun::Parameterisation> > strfuns;
  vector<TGraph*> g_strfuns_f2, g_strfuns_fl;
  for (const auto& sf_type : strfun_types) {
    auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(sf_type);
    std::ostringstream oss;
    oss << (cepgen::strfun::Type)sf->name();
    const auto sf_name = oss.str();
    g_strfuns_f2.emplace_back(new TGraph);
    (*g_strfuns_f2.rbegin())->SetTitle((sf_name + ";x_{Bj};F_{2}").c_str());
    g_strfuns_fl.emplace_back(new TGraph);
    (*g_strfuns_fl.rbegin())->SetTitle((sf_name + ";x_{Bj};F_{L}").c_str());
    strfuns.emplace_back(move(sf));
  }
  for (int i = 0; i < num_points; ++i) {
    const double x = xmin + i * (xmax - xmin) / (num_points - 1);
    out << x << "\t";
    size_t j = 0;
    for (auto& sf : strfuns) {
      const auto strfun = (*sf)(x, q2);
      out << "\t" << strfun.F2 << "\t" << strfun.FL;
      g_strfuns_f2.at(j)->SetPoint(g_strfuns_f2.at(j)->GetN(), x, strfun.F2);
      g_strfuns_fl.at(j)->SetPoint(g_strfuns_fl.at(j)->GetN(), x, strfun.FL);
      ++j;
    }
    out << "\n";
  }
  CG_LOG << "Scan written in \"" << output_file << "\".";
  out.close();

  for (auto& plt : map<const char*, vector<TGraph*> >{{"F2", g_strfuns_f2}, {"FL", g_strfuns_fl}}) {
    cepgen::ROOTCanvas c(plt.first, Form("Q^{2} = %g GeV^{2}", q2));
    c.SetLogy();
    TMultiGraph mg;
    size_t i = 0;
    for (auto& gr : plt.second) {
      mg.Add(gr);
      gr->SetLineColor(cepgen::ROOTCanvas::colours[i]);
      c.AddLegendEntry(gr, gr->GetTitle(), "l");
      ++i;
    }
    mg.Draw("al");
    // ugly fix to propagate first plot axes label onto multigraph
    mg.GetHistogram()->SetTitle((*plt.second.begin())->GetTitle());
    c.Prettify(mg.GetHistogram());
    c.Save("pdf");
  }

  return 0;
}
