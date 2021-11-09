/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#include <TGraph2D.h>

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_card;
  int npoints;
  vector<int> dim;
  double def;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("input,i", "input card", &input_card)
      .addOptionalArgument("default,D", "default value for non-varying coordinates", &def, 0.5)
      .addOptionalArgument("dim,d", "dimensions to probe", &dim, vector<int>{})
      .addOptionalArgument("num-points,n", "number of points to probe", &npoints, 100)
      .parse();

  TGraph gr_scan_1d;
  TGraph2D gr_scan_2d;
  if (dim.size() > 3)
    throw CG_FATAL("main") << "Number of dimensions to probe (" << dim.size() << ") is too high";

  cepgen::Generator gen;
  gen.setParameters(cepgen::card::Handler::parse(input_card));
  CG_LOG << gen.parameters();
  const size_t ndim = gen.parameters()->process().ndim();

  vector<double> coord(ndim, def);

  for (int i = 0; i < npoints; ++i) {
    const double x = i * 1. / npoints;
    switch (dim.size()) {
      case 0:
        gr_scan_1d.SetPoint(gr_scan_1d.GetN(), x, gen.computePoint(vector<double>(ndim, x)));
        break;
      case 1:
        coord[dim.at(0)] = x;
        gr_scan_1d.SetPoint(gr_scan_1d.GetN(), x, gen.computePoint(coord));
        break;
      case 2:
        coord[dim.at(0)] = x;
        for (int j = 0; j < npoints; ++j) {
          const double y = j * 1. / npoints;
          coord[dim.at(1)] = y;
          gr_scan_2d.SetPoint(gr_scan_2d.GetN(), x, y, gen.computePoint(coord));
        }
        break;
    }
  }
  {
    cepgen::ROOTCanvas c("test_scan");
    gStyle->SetPalette(kBeach);
    string xlabel, ylabel;
    switch (dim.size()) {
      case 0:
      case 1: {
        if (dim.empty())
          xlabel = Form("x_{i = 0, ..., %ld}", ndim - 1);
        else
          xlabel = Form("x_{%d}", dim.at(0));
        gr_scan_1d.SetMarkerStyle(24);
        c.SetTopLabel(Form("%s variation, all others x_{i} at %g", xlabel.c_str(), def));
        gr_scan_1d.SetTitle(Form(";%s;d^{%ld}#sigma/d#bf{x}^{%ld}", xlabel.c_str(), ndim, ndim));
        gr_scan_1d.Draw("ap");
        c.Prettify(gr_scan_1d.GetHistogram());
        c.SetLogy();
      } break;
      case 2: {
        xlabel = Form("x_{%d}", dim.at(0));
        ylabel = Form("x_{%d}", dim.at(1));
        c.SetTopLabel(Form("(%s, %s) variation, all others x_{i} at %g", xlabel.c_str(), ylabel.c_str(), def));
        gr_scan_2d.SetTitle(Form(";%s;%s;d^{%ld}#sigma/d#bf{x}^{%ld}", xlabel.c_str(), ylabel.c_str(), ndim, ndim));
        gr_scan_2d.Draw("colz");
        c.Prettify((TH1*)gr_scan_2d.GetHistogram());
        c.SetLogz();
      }
    }
    c.Save("pdf");
  }

  return 0;
}
