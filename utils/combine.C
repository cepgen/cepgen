//#include "style.C"

void combine()
{
  //SetStyle();

  ifstream in;

  const Int_t lepPdg = 13;
  const Int_t n=100;
  const Int_t maxEvts = -1;
  //const Int_t maxEvts = 5e5;
  const Double_t ep = 3500.;
  const Int_t nHist = 14;

  stringstream ss;
  Int_t i, j;

  Double_t px, py, pz, pt, e, m, eta, q2m, pp3, pp5, weight;
  Int_t pdg;

  TLine *line;
  TH1D *h[nHist], *h_2[nHist];
  TH1D *htmp;
  TCanvas *c[nHist];
  TLegend *leg;
  TLorentzVector lep1, lep2, prot;
  
  // LPAIR's TTree definition
  Float_t px_[n], py_[n], pz_[n], e_[n], m_[n], eta_[n];
  Int_t pdgId_[n], npart_;
  TFile *lp;
  TTree *tree;
  bool lep1set, lep2set, pset;

  //lp = new TFile("/home/forthomme/LpairAnalysis/trunk/samples/lpair-pt5-mumu-elastic.root");
  lp = new TFile("/home/forthomme/LpairAnalysis/trunk/samples/lpair-7tev-elastic-pt5.root");
  //lp = new TFile("/home/forthomme/LpairAnalysis/trunk/samples/lpair-mumu-pt15-8tev-elastic.root");
  //lp = new TFile("/home/forthomme/LpairAnalysis/trunk/samples/lpair-7tev-elastic-nocuts.root");
  //lp = new TFile("/home/forthomme/LpairAnalysis/trunk/samples/lpair-7tev-elastic-pt10-theta0to180.root");
  tree = (TTree*)(lp->Get("h4444"));
  tree->SetBranchAddress("px", px_);
  tree->SetBranchAddress("py", py_);
  tree->SetBranchAddress("pz", pz_);
  tree->SetBranchAddress("E", e_);
  tree->SetBranchAddress("m", m_);
  tree->SetBranchAddress("Eta", eta_);
  tree->SetBranchAddress("icode", pdgId_);
  tree->SetBranchAddress("ip", &npart_);

  gStyle->SetOptStat(0);

  h[0] = new TH1D("pt", "p_{T}", 200, 0., 100.);
  h_2[0] = new TH1D("pt_2", "p_{T}", 200, 0., 100.);
  h[1] = new TH1D("px", "p_{x}", 500, -100., 100.);
  h_2[1] = new TH1D("px_2", "p_{x}", 500, -100., 100.);
  h[2] = new TH1D("py", "p_{y}", 500, -100., 100.);
  h_2[2] = new TH1D("py_2", "p_{y}", 500, -100., 100.);
  h[3] = new TH1D("pz", "p_{z}", 200, -100., 100.);
  h_2[3] = new TH1D("pz_2", "p_{z}", 200, -100., 100.);
  h[4] = new TH1D("e", "E", 200, 0., 100.);
  h_2[4] = new TH1D("e_2", "E", 200, 0., 100.);
  h[5] = new TH1D("p", "p", 200, 0., 100.);
  h_2[5] = new TH1D("p_2", "p", 200, 0., 100.);
  h[6] = new TH1D("eta", "#eta", 200, -10., 10.);
  h_2[6] = new TH1D("eta_2", "#eta", 200, -10., 10.);
  h[7] = new TH1D("phi", "#phi", 200, -5., 5.);
  h_2[7] = new TH1D("phi_2", "#phi", 200, -5., 5.);
  h[8] = new TH1D("theta", "#theta", 140, 0., 3.5);
  h_2[8] = new TH1D("theta_2", "#theta", 140, 0., 3.5);
  h[9] = new TH1D("mass", "m(l^{+}l^{-})", 200, 0., 100.);
  h_2[9] = new TH1D("mass_2", "m(l^{+}l^{-})", 200, 0., 100.);
  h[10] = new TH1D("ptpair", "p_{T}(l^{+}l^{-})", 100, 0., 5.);
  h_2[10] = new TH1D("ptpair_2", "p_{T}(l^{+}l^{-})", 100, 0., 5.);
  h[11] = new TH1D("ypair", "y(l^{+}l^{-})", 100, -15., 15.);
  h_2[11] = new TH1D("ypair_2", "y(l^{+}l^{-})", 100, -15., 15.);
  h[12] = new TH1D("q2m", "Q^{2}", 200, 0., 100.);
  h_2[12] = new TH1D("q2m_2", "Q^{2}", 200, 0., 100.);
  h[13] = new TH1D("pp", "p_{proton}", (Int_t)ep/20, 0., ep);
  h_2[13] = new TH1D("pp_2", "p_{proton}", (Int_t)ep/20, 0., ep);
  /*h[11] = new TH1D("pp5", "p_{5}", 100, 0., 10.);
  h_2[11] = new TH1D("pp5_2", "p_{5}", 100, 0., 10.);*/

  // First fetch the LPAIR++ output
  in.open("test");
  i = 0;
  lep1set = lep2set = false;
  while(in >> e >> px >> py >> pz >> pt >> m >> eta >> pdg >> weight) {
    if (i%2==0 && (i/2)%10000==0) {
      cout << "[LPAIR++] Event #" << i/2 << endl;
    }
    if (i<5)
      cout << i 
	   << "\t" << pdg
	   << "\t" << m
	   << "\t" << eta
	   << "\t" << px
	   << "\t" << py
	   << "\t" << pz
	   << "\t" << pt
	   << "\t" << e
	   << endl;
    if (pdg>0) {
      lep1.SetXYZM(px, py, pz, m);
      lep1set = true;
    }
    else {
      lep2.SetXYZM(px, py, pz, m);
      lep2set = true;
    }
    if (lep1set && lep2set) {
      //if (fabs((lep1+lep2).Rapidity())>5.) continue;
      h[0]->Fill(lep1.Pt());
      h[1]->Fill(lep1.Px());
      h[2]->Fill(lep1.Py());
      h[3]->Fill(lep1.Pz());
      h[4]->Fill(lep1.E());
      h[5]->Fill(lep1.P());
      h[6]->Fill(eta);
      h[7]->Fill(lep1.Phi());
      h[8]->Fill(lep1.Theta());
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
      h[9]->Fill((lep1+lep2).M());
      h[10]->Fill((lep1+lep2).Pt());
      h[11]->Fill((lep1+lep2).Rapidity());
      lep1set = lep2set = false;
    }
    if (maxEvts>0 && i>maxEvts) break;
    i++;
  }
  in.close();

  in.open("test_q2");
  while (in >> q2m >> pp3 >> pp5) {
    h[12]->Fill(-q2m);
    h[13]->Fill(pp3);
    h[13]->Fill(pp5);
  }
  in.close();


  // Then fetch the LPAIR output (converted as a TTree)
  for (i=0; i<tree->GetEntries(); i++) {
    if (i%10000==0) {
      cout << "[ LPAIR ] Event #" << i << endl;
    }
    lep1set = lep2set = false;
    pset = false;
    tree->GetEntry(i);
    for (j=0; j<npart_; j++) {
      if (abs(pdgId_[j])==2212) {
        prot.SetXYZM(px_[j], py_[j], pz_[j], m_[j]);
        h_2[13]->Fill(prot.P());
      }
      if (abs(pdgId_[j])==2212 && !pset) {
        prot.SetXYZM(px_[j], py_[j], pz_[j], m_[j]);
        q2m = -(prot.P()-ep);
        h_2[12]->Fill(q2m);
        pset = true;
      }
      if (abs(pdgId_[j])!=lepPdg) continue;
      if (pdgId_[j]>0) {
	lep1.SetXYZM(px_[j], py_[j], pz_[j], m_[j]);
	lep1set = true;
      }
      else {
	lep2.SetXYZM(px_[j], py_[j], pz_[j], m_[j]);
	lep2set = true;
      }
    }
    if (lep1set && lep2set) {
      //if (fabs((lep1+lep2).Rapidity())>5.) continue;
      h_2[0]->Fill(lep1.Pt());
      h_2[1]->Fill(lep1.Px());
      h_2[2]->Fill(lep1.Py());
      h_2[3]->Fill(lep1.Pz());
      h_2[4]->Fill(lep1.E());
      h_2[5]->Fill(lep1.P());
      h_2[6]->Fill(lep1.Eta());
      h_2[7]->Fill(lep1.Phi());
      h_2[8]->Fill(lep1.Theta());
      h_2[9]->Fill((lep1+lep2).M());
      h_2[10]->Fill((lep1+lep2).Pt());
      h_2[11]->Fill((lep1+lep2).Rapidity());
      lep1set = lep2set = false;
    }
    if (maxEvts>0 && i>maxEvts) break;
  }

  leg = new TLegend(.82, .65, .95, .75);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kBlack);

  for (i=0; i<nHist; i++) {
    c[i] = new TCanvas(); 
    c[i]->Divide(1,2);

    TPad *c_1 = (TPad*)(c[i]->GetPad(1));
    c_1->SetPad(0.,.250,1.,1.);
    c_1->SetRightMargin(0.03);
    c_1->SetBottomMargin(0.);
    c_1->SetGrid(1,1);
    TPad *c_2 = (TPad*)(c[i]->GetPad(2));
    c_2->SetPad(0.,0.,1.,.250);
    c_2->SetBottomMargin(0.3);
    c_2->SetRightMargin(0.03);
    c_2->SetTopMargin(0.);

    c[i]->cd(1);
    h[i]->SetFillColor(kBlue);
    h[i]->SetFillStyle(3005);
    h[i]->SetLineColor(kBlack);
    //ss.str(""); ss << "#frac{1}{#sigma} #frac{d#sigma}{d" << h[i]->GetTitle() << "}";
    ss.str(""); ss << "#frac{dN}{d" << h[i]->GetTitle() << "}";
    h_2[i]->GetYaxis()->SetTitle(ss.str().c_str());
    h[i]->GetXaxis()->SetTitleFont(43);
    h[i]->GetXaxis()->SetTitleSize(14);
    h[i]->GetXaxis()->SetTitleOffset(4.);
    h[i]->GetYaxis()->SetTitleFont(43);
    h[i]->GetYaxis()->SetTitleSize(14);
    h[i]->GetYaxis()->SetTitleOffset(1.4);
    h_2[i]->GetYaxis()->SetTitleFont(43);
    h_2[i]->GetYaxis()->SetTitleSize(14);
    h_2[i]->GetYaxis()->SetTitleOffset(1.4);
    h[i]->GetXaxis()->SetLabelFont(43);
    h[i]->GetXaxis()->SetLabelSize(14);
    h[i]->GetYaxis()->SetLabelFont(43);
    h[i]->GetYaxis()->SetLabelSize(14);
    //h[i]->Scale(1./h[i]->Integral());
    h_2[i]->SetFillColor(kRed);
    h_2[i]->SetFillStyle(3004);
    h_2[i]->SetLineColor(kBlack);
    //h_2[i]->Scale(1./h_2[i]->Integral());
    h_2[i]->Draw();
    h[i]->Draw("SAME");
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
    htmp->GetXaxis()->SetTitle(h[i]->GetTitle());
    htmp->GetYaxis()->SetTitle("LPAIR++/LPAIR");
    htmp->Draw("E");
    line->Draw();
  }
}
