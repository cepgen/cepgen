// clang-format off
R__ADD_INCLUDE_PATH(/home/forthomme/work/dev/cepgen)
R__ADD_LIBRARY_PATH(../../../build)
R__LOAD_LIBRARY(libCepGen.so)
R__LOAD_LIBRARY(libCepGenRoot.so)
// clang-format on

#include "CepGen/Generator.h"
#include "CepGenAddOns/ROOTWrapper/utils/hist_utils.h"

using namespace std;

void compare_root_trees(const string& base,
                        const string& comp,
                        const string& base_label = "",
                        const string& comp_label = "") {
  cepgen::initialise();
  auto hists_base = fill_histograms(base);
  auto hists_comp = fill_histograms(comp);

  for (size_t i = 0; i < hists_base.size(); ++i) {
    cepgen::ROOTCanvas c(hists_base[i]->GetName());
    THStack hs;
    hists_base[i]->SetLineColor(cepgen::ROOTCanvas::colours[0]);
    hists_base[i]->Scale(1. / hists_comp[i]->Integral());
    if (!base_label.empty())
      c.AddLegendEntry(hists_base[i], base_label.c_str(), "l");
    hists_comp[i]->SetLineColor(cepgen::ROOTCanvas::colours[1]);
    hists_comp[i]->Scale(1. / hists_comp[i]->Integral());
    if (!comp_label.empty())
      c.AddLegendEntry(hists_comp[i], comp_label.c_str(), "l");
    hs.Add(hists_base[i]);
    hs.Add(hists_comp[i]);
    hs.Draw("hist,nostack");
    hs.GetHistogram()->GetXaxis()->SetTitle(hists_base[i]->GetXaxis()->GetTitle());
    hs.GetHistogram()->GetYaxis()->SetTitle(hists_base[i]->GetYaxis()->GetTitle());
    c.Prettify(hs.GetHistogram());
    c.Save("pdf");
  }
}
