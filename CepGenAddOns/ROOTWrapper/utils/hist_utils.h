#ifndef CepGenAddOns_ROOTWrapper_utils_hist_utils_h
#define CepGenAddOns_ROOTWrapper_utils_hist_utils_h

#include "CepGen/Event/Event.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"
#include "CepGenAddOns/ROOTWrapper/ROOTTreeInfo.h"

typedef vector<TH1D*> hists_t;

hists_t fill_histograms(const string& filename) {
  hists_t out = {new TH1D("invm", ";m_{central} (GeV);d#sigma/dm", 200, 150., 550.),
                 new TH1D("ptpair", ";p_{T}^{central} (GeV);d#sigma/dp_{T}", 100, 0., 5.),
                 new TH1D("acop", ";1-|#Delta#phi/#pi|;d#sigma/d#Delta#phi)", 50, 0., 1.e-2)};
  auto file = TFile::Open(filename.c_str(), "r");
  ROOT::CepGenRun run;
  run.attach(file);
  cout << ">>> " << run.process_name << ": " << run.process_parameters << endl;

  ROOT::CepGenEvent evt_tree;
  evt_tree.attach(file);
  cepgen::Event evt;
  while (evt_tree.next(evt)) {
    out[0]->Fill(evt[4].mass());
    out[1]->Fill(evt[4].momentum().pt());
    //out[2]->Fill(
  }
  return out;
}

#endif
