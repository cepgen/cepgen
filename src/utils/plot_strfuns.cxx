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
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<int> strfun_types;
  double q2, xmin, xmax;
  int var, num_points;
  string output_file;
  bool logx, logy;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("sf,s", "structure functions modelling", &strfun_types)
      .addOptionalArgument("q2,q", "parton virtuality (GeV^2)", &q2, 10.)
      .addOptionalArgument("var,t", "variable to study (0=xBj, 1=w)", &var, 0)
      .addOptionalArgument("xmax,m", "minimal Bjorken x", &xmin, 1.e-7)
      .addOptionalArgument("xmax,M", "maximal Bjorken x", &xmax, 1.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "strfuns.scan.output.txt")
      .addOptionalArgument("logx", "logarithmic x-axis", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-axis", &logy, false)
      .parse();

  const double lxmin = log10(xmin), lxmax = log10(xmax);

  cepgen::initialise();

  string var_name;
  switch (var) {
    case 0:
      var_name = "x_{Bj}";
      break;
    case 1:
      var_name = "w (GeV)";
      break;
    case 2:
      var_name = "w^{2} (GeV^{2})";
      break;
    default:
      throw CG_FATAL("main") << "Unsupported variable to be plotted!";
  }

  ofstream out(output_file);
  out << "# structure functions: ";
  string sep;
  for (const auto& sf_type : strfun_types)
    out << sep << sf_type, sep = ", ";
  out << "\n"
      << "# x in [" << xmin << ", " << xmax << "]\n";

  const float mp = cepgen::PDG::get().mass(2212), mp2 = mp * mp;

  vector<unique_ptr<cepgen::strfun::Parameterisation> > strfuns;
  vector<TGraph*> g_strfuns_f2, g_strfuns_fl, g_strfuns_fe, g_strfuns_fm, g_strfuns_w1, g_strfuns_w2;
  for (const auto& sf_type : strfun_types) {
    auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(sf_type);
    const auto sf_name = cepgen::strfun::StructureFunctionsFactory::get().describe(sf_type);
    g_strfuns_f2.emplace_back(new TGraph);
    (*g_strfuns_f2.rbegin())->SetTitle((sf_name + ";" + var_name + ";F_{2}(" + var_name + ", Q^{2})").c_str());
    g_strfuns_fl.emplace_back(new TGraph);
    (*g_strfuns_fl.rbegin())->SetTitle((sf_name + ";" + var_name + ";F_{L}(" + var_name + ", Q^{2})").c_str());
    g_strfuns_fe.emplace_back(new TGraph);
    (*g_strfuns_fe.rbegin())->SetTitle((sf_name + ";" + var_name + ";F_{E}(" + var_name + ", Q^{2})").c_str());
    g_strfuns_fm.emplace_back(new TGraph);
    (*g_strfuns_fm.rbegin())->SetTitle((sf_name + ";" + var_name + ";F_{M}(" + var_name + ", Q^{2})").c_str());
    g_strfuns_w1.emplace_back(new TGraph);
    (*g_strfuns_w1.rbegin())->SetTitle((sf_name + ";" + var_name + ";W_{1}(" + var_name + ", Q^{2})").c_str());
    g_strfuns_w2.emplace_back(new TGraph);
    (*g_strfuns_w2.rbegin())->SetTitle((sf_name + ";" + var_name + ";W_{2}(" + var_name + ", Q^{2})").c_str());
    strfuns.emplace_back(move(sf));
  }
  for (int i = 0; i < num_points; ++i) {
    const double x =
        (!logx) ? xmin + i * (xmax - xmin) / (num_points - 1) : pow(10, lxmin + i * (lxmax - lxmin) / (num_points - 1));
    out << x << "\t";
    size_t j = 0;
    for (auto& sf : strfuns) {
      double xbj;
      switch (var) {
        case 0:
          xbj = x;
          break;
        case 1:
          xbj = cepgen::utils::xBj(q2, mp2, x * x);
          break;
        case 2:
          xbj = cepgen::utils::xBj(q2, mp2, x);
          break;
        default:
          xbj = 0.;
          break;
      }
      out << "\t" << sf->F2(xbj, q2) << "\t" << sf->FL(xbj, q2);
      g_strfuns_f2.at(j)->SetPoint(g_strfuns_f2.at(j)->GetN(), x, sf->F2(xbj, q2));
      g_strfuns_fl.at(j)->SetPoint(g_strfuns_fl.at(j)->GetN(), x, sf->FL(xbj, q2));
      g_strfuns_fe.at(j)->SetPoint(g_strfuns_fe.at(j)->GetN(), x, sf->FE(xbj, q2));
      g_strfuns_fm.at(j)->SetPoint(g_strfuns_fm.at(j)->GetN(), x, sf->FM(xbj, q2));
      g_strfuns_w1.at(j)->SetPoint(g_strfuns_w1.at(j)->GetN(), x, sf->W1(xbj, q2));
      g_strfuns_w2.at(j)->SetPoint(g_strfuns_w2.at(j)->GetN(), x, sf->W2(xbj, q2));
      ++j;
    }
    out << "\n";
  }
  CG_LOG << "Scan written in \"" << output_file << "\".";
  out.close();

  for (auto& plt : map<string, vector<TGraph*> >{{"f2", g_strfuns_f2},
                                                 {"fl", g_strfuns_fl},
                                                 {"fe", g_strfuns_fe},
                                                 {"fm", g_strfuns_fm},
                                                 {"w1", g_strfuns_w1},
                                                 {"w2", g_strfuns_w2}}) {
    cepgen::ROOTCanvas c(("sfcomp_" + plt.first).c_str(), Form("Q^{2} = %g GeV^{2}", q2));
    TMultiGraph mg;
    if (logx)
      c.SetLogx();
    if (logy)
      c.SetLogy();
    size_t i = 0;
    for (auto& gr : plt.second) {
      mg.Add(gr);
      gr->SetLineColor(cepgen::ROOTCanvas::colours[i]);
      c.AddLegendEntry(gr, gr->GetTitle(), "l");
      ++i;
    }
    mg.Draw("al");
    // ugly fix to propagate first plot axes label onto multigraph
    mg.GetHistogram()->GetXaxis()->SetTitle((*plt.second.begin())->GetXaxis()->GetTitle());
    mg.GetHistogram()->GetYaxis()->SetTitle((*plt.second.begin())->GetYaxis()->GetTitle());
    c.Prettify(mg.GetHistogram());
    c.Save("pdf");
  }

  return 0;
}
