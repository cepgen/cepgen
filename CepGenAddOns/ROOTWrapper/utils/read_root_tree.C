// clang-format off
R__ADD_INCLUDE_PATH(/home/forthomme/work/dev/cepgen)
R__ADD_LIBRARY_PATH(../../../build)
R__LOAD_LIBRARY(libCepGen.so)
R__LOAD_LIBRARY(libCepGenRoot.so)
// clang-format on

#include "CepGen/Generator.h"
#include "CepGenAddOns/ROOTWrapper/utils/hist_utils.h"

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
