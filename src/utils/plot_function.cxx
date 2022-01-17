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

#include <TF1.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Functional.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

using namespace std;

int main(int argc, char* argv[]) {
  string function;
  int num_points;
  double min_x, max_x;
  bool draw_grid;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("function,f", "function to parse", &function, "min(1,exp(-x/10))")
      .addOptionalArgument("num-points,n", "number of points to consider", &num_points, 100)
      .addOptionalArgument("min-x,l", "minimal range", &min_x, -5.)
      .addOptionalArgument("max-x,H", "maximal range", &max_x, +5.)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .parse();

  TGraph gr_rt;
  TF1 f_rt("f_rt", "TMath::Min(1.,TMath::Exp(-x/10))", min_x, max_x);
  for (unsigned short i = 0; i < num_points; ++i) {
    const double x = min_x + (max_x - min_x) / (num_points - 1) * i;
    gr_rt.SetPoint(i, x, f_rt.Eval(x));
  }

  CG_LOG << "Function to be plotted: " << function;

  map<string, TGraph> m_gr_fb, m_gr_diff;
  for (const auto& func : cepgen::utils::FunctionalFactory::get().modules()) {
    CG_LOG << "Building \"" << func << "\" functional.";
    try {
      auto test =
          cepgen::utils::FunctionalFactory::get().build(func,
                                                        cepgen::ParametersList()
                                                            .set<std::string>("expression", function)
                                                            .set<std::vector<std::string> >("variables", {"x"}));
      double chi2 = 0;
      for (unsigned short i = 0; i < num_points; ++i) {
        const double x = min_x + (max_x - min_x) / (num_points - 1) * i;
        const double val = (*test)(x), val_ref = gr_rt.GetY()[i];
        m_gr_fb[func].SetPoint(i, x, val);
        m_gr_diff[func].SetPoint(i, x, val - val_ref);
        chi2 += pow(val - val_ref, 2);
      }
      chi2 = sqrt(chi2);
      if (chi2 > 1.e-9)
        throw CG_FATAL("main") << "Test failed with chi2 = " << chi2 << "!";
    } catch (const cepgen::Exception&) {
      CG_WARNING("main") << "Exception encountered in \"" << func << "\" functional builder.";
      continue;
    }
  }

  CG_LOG << "Test passed!";

  {
    cepgen::ROOTCanvas c("test_graph", "CepGen validation", true);
    if (draw_grid)
      c.SetGrid(true, true);
    TMultiGraph mg;
    mg.Add(&gr_rt);
    c.AddLegendEntry(&gr_rt, "ROOT", "l");
    size_t i = 0;
    for (auto& gr_fb : m_gr_fb) {
      mg.Add(&gr_fb.second);
      gr_fb.second.SetLineWidth(3);
      gr_fb.second.SetLineStyle(2 + (i++));
      c.AddLegendEntry(&gr_fb.second, Form("Functional (%s)", gr_fb.first.c_str()), "l");
    }
    i = 0;
    for (auto& gr_diff : m_gr_diff) {
      gr_diff.second.SetLineStyle(2 + i);
      gr_diff.second.SetLineColor(cepgen::ROOTCanvas::colours[i]);
      gr_diff.second.Draw("same");
      ++i;
    }
    mg.Draw("al");
    c.Prettify(mg.GetHistogram());
    c.Save("pdf");
    i = 0;
    for (auto& gr_fb : m_gr_fb) {
      auto* ratio =
          c.RatioPlot(gr_fb.second.GetHistogram(), {gr_rt.GetHistogram()}, -1., 1., (i == 0 ? "al" : "l,same"))[0];
      ratio->SetLineColor(kRed);
      ratio->SetLineWidth(3);
      ratio->SetLineStyle(2 + (i++));
    }
  }
  return 0;
}
