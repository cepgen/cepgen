/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#ifndef CepGenRoot_utils_hist_utils_h
#define CepGenRoot_utils_hist_utils_h

#include "CepGen/Event/Event.h"
#include "CepGenRoot/ROOTCanvas.h"
#include "CepGenRoot/ROOTTreeInfo.h"

typedef vector<TH1D*> hists_t;

hists_t fill_histograms(const string& filename) {
  hists_t out = {new TH1D("invm", ";m_{central} (GeV);d#sigma/dm", 200, 0., 20.),
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
    evt.dump();
    out[0]->Fill(evt(4).momentum().mass());
    out[1]->Fill(evt(4).momentum().pt());
    out[2]->Fill(1. - fabs(evt(7).momentum().deltaPhi(evt(8).momentum()) * M_1_PI));
  }
  return out;
}

#endif
