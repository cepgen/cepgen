/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

// clang-format off
R__ADD_INCLUDE_PATH(/home/forthomme/work/dev/cepgen)
R__ADD_LIBRARY_PATH(../../../build)
R__LOAD_LIBRARY(libCepGen.so)
R__LOAD_LIBRARY(libCepGenRoot.so)
// clang-format on

#include "CepGen/Generator.h"
#include "hist_utils.h"

using namespace std;

void read_root_tree(const string& base, const string& label = "") {
  cepgen::initialise();
  auto hists = fill_histograms(base);

  for (size_t i = 0; i < hists.size(); ++i) {
    cepgen::ROOTCanvas c(hists[i]->GetName());
    THStack hs;
    hists[i]->SetLineColor(cepgen::ROOTCanvas::colours[0]);
    if (!label.empty())
      c.AddLegendEntry(hists[i], label.c_str(), "l");
    hists[i]->Draw("hist");
    c.Prettify(hists[i]);
    c.Save("pdf");
  }
}
