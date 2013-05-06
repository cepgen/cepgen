void combine()
{
  ifstream in;
  in.open("test");

  const Int_t lepPdg = 13;
  const Int_t n=100;

  stringstream ss;
  Int_t i, j;

  Double_t px, py, pz, pt, e, m;
  Int_t pdg;

  TLine *line;
  TH1D *h[7], *h_2[7];
  TH1D *htmp;
  TCanvas *c[7];
  TLegend *leg;
  TLorentzVector lep1, lep2;
  
  // LPAIR's TTree definition
  Float_t px_[n], py_[n], pz_[n], e_[n], m_[n];
  Int_t pdgId_[n], npart_;
  TFile *lp;
  TTree *tree;
  bool lep1set, lep2set;

  lp = new TFile("lpair-pt5-mumu-elastic.root");
  tree = (TTree*)(lp->Get("h4444"));
  tree->SetBranchAddress("px", px_);
  tree->SetBranchAddress("py", py_);
  tree->SetBranchAddress("pz", pz_);
  tree->SetBranchAddress("E", e_);
  tree->SetBranchAddress("m", m_);
  tree->SetBranchAddress("icode", pdgId_);
  tree->SetBranchAddress("ip", &npart_);

  h[0] = new TH1D("pt", "p_{T}", 200, 0., 100.);
  h_2[0] = new TH1D("pt_2", "p_{T}", 200, 0., 100.);
  h[1] = new TH1D("px", "p_{x}", 200, -100., 100.);
  h_2[1] = new TH1D("px_2", "p_{x}", 200, -100., 100.);
  h[2] = new TH1D("py", "p_{y}", 200, -100., 100.);
  h_2[2] = new TH1D("py_2", "p_{y}", 200, -100., 100.);
  h[3] = new TH1D("pz", "p_{z}", 200, -100., 100.);
  h_2[3] = new TH1D("pz_2", "p_{z}", 200, -100., 100.);
  h[4] = new TH1D("e", "E", 200, 0., 100.);
  h_2[4] = new TH1D("e_2", "E", 200, 0., 100.);
  h[5] = new TH1D("mass", "m(l^{+}l^{-})", 200, 0., 100.);
  h_2[5] = new TH1D("mass_2", "m(l^{+}l^{-})", 200, 0., 100.);
  h[6] = new TH1D("ptpair", "p_{T}(l^{+}l^{-})", 100, 0., 5.);
  h_2[6] = new TH1D("ptpair_2", "p_{T}(l^{+}l^{-})", 100, 0., 5.);

  // First fetch the LPAIR++ output
  i = 0;
  lep1set = lep2set = false;
  while(!in.eof()) {
    in >> e >> px >> py >> pz >> pt >> m >> pdg;
    if (i%2==0 && (i/2)%10000==0) {
      cout << "[LPAIR++] Event #" << i/2 << endl;
    }
    if (pdg>0) {
      //lep1.SetPxPyPzE(px, py, pz, e);
      lep1.SetXYZM(px, py, pz, m);
      h[0]->Fill(pt);
      h[1]->Fill(px);
      h[2]->Fill(py);
      h[3]->Fill(pz);
      h[4]->Fill(e);
      lep1set = true;
    }
    else {
      //lep2.SetPxPyPzE(px, py, pz, e);
      lep2.SetXYZM(px, py, pz, m);
      lep2set = true;
    }
    if (lep1set && lep2set) {
      /*double mass = std::sqrt(
			      std::pow(lep1.M(),2)+
			      std::pow(lep2.M(),2)+
			      2*(
				 lep1.E()*lep2.E()-
				 lep1.Px()*lep2.Px()-
				 lep1.Py()*lep2.Py()-
				 lep1.Pz()*lep2.Pz()
				 )
			      );
			      hmass->Fill(mass);*/
      h[5]->Fill((lep1+lep2).M());
      h[6]->Fill((lep1+lep2).Pt());
      lep1set = lep2set = false;
    }
    if (i>2e6) break;
    i++;
  }

  // Then fetch the LPAIR output (converted as a TTree)
  for (i=0; i<tree->GetEntries(); i++) {
    if (i%10000==0) {
      cout << "[ LPAIR ] Event #" << i << endl;
    }
    lep1set = lep2set = false;
    tree->GetEntry(i);
    for (j=0; j<npart_; j++) {
      if (abs(pdgId_[j])!=lepPdg) continue;
      if (pdgId_[j]>0) {
	lep1.SetXYZM(px_[j], py_[j], pz_[j], m_[j]);
	h_2[0]->Fill(lep1.Pt());
	h_2[1]->Fill(px_[j]);
	h_2[2]->Fill(py_[j]);
	h_2[3]->Fill(pz_[j]);
	h_2[4]->Fill(e_[j]);
	lep1set = true;
      }
      else {
	lep2.SetXYZM(px_[j], py_[j], pz_[j], m_[j]);
	lep2set = true;
      }
    }
    if (lep1set && lep2set) {
      h_2[5]->Fill((lep1+lep2).M());
      h_2[6]->Fill((lep1+lep2).Pt());
    }
  }

  leg = new TLegend(.82, .65, .95, .75);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kBlack);

  for (i=0; i<7; i++) {
    c[i] = new TCanvas(); 
    c[i]->Divide(1,2);
    c[i]->cd(1);
    h[i]->SetFillColor(kBlue);
    h[i]->SetFillStyle(3005);
    h[i]->SetLineColor(kBlack);
    ss.str(""); ss << "#frac{1}{#sigma} #frac{d#sigma}{d" << h[i]->GetTitle() << "}";
    h[i]->GetYaxis()->SetTitle(ss.str().c_str());
    h[i]->Scale(1./h[i]->Integral());
    h[i]->Draw();
    h_2[i]->SetFillColor(kRed);
    h_2[i]->SetFillStyle(3004);
    h_2[i]->SetLineColor(kBlack);
    h_2[i]->Scale(1./h_2[i]->Integral());
    h_2[i]->Draw("SAME");
    if (i==0) {
      leg->AddEntry(h[i], "LPAIR++");
      leg->AddEntry(h_2[i], "LPAIR");
    }
    leg->Draw("SAME");
    c[i]->cd(2);
    line = new TLine(h[i]->GetXaxis()->GetXmin(),1., h[i]->GetXaxis()->GetXmax(), 1.);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    htmp = (TH1D*)h[i]->Clone();
    htmp->Divide(h_2[i]);
    htmp->SetTitle("");
    htmp->GetYaxis()->SetTitle("LPAIR++/LPAIR");
    htmp->Draw("E");
    line->Draw();
  }
}
