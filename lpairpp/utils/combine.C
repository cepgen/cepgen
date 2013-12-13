//#include "style.C"

#define PTSINGLE 0
#define PXSINGLE 1
#define PYSINGLE 2
#define PZSINGLE 3
#define ESINGLE 4
#define PSINGLE 5
#define ETASINGLE 6
#define PHISINGLE 7
#define THETASINGLE 8
#define DPT 9
#define ACOP 10
#define MPAIR 11
#define PTPAIR 12
#define YPAIR 13
#define Q2 14
#define PPROTON 15
#define T1 16
#define T1MIN 17
#define T1MAX 18
#define T2 19
#define S1 20
#define S2 21
#define D3 22
#define WTREAT 23
#define ZTREAT 24
#define XIN0TREAT 25
#define XIN1TREAT 26
#define XIN2TREAT 27
#define XIN3TREAT 28
#define XIN4TREAT 29
#define XIN5TREAT 30
#define XIN6TREAT 31
#define XOUT0TREAT 32
#define XOUT1TREAT 33
#define XOUT2TREAT 34
#define XOUT3TREAT 35
#define XOUT4TREAT 36
#define XOUT5TREAT 37
#define XOUT6TREAT 38

#define NHIST 39

void combine()
{
  //SetStyle();

  ifstream in;

  const Int_t lepPdg = 13;
  const Int_t n=100;
  //const Int_t maxEvts = -1;
  const Int_t maxEvts = 2e4;
  const Double_t ep = 3500.;
  const Double_t pi = acos(-1);

  stringstream ss;
  Int_t i, j;

  Double_t px, py, pz, pt, e, m, eta, q2m, pp3, pp5, weight;
  Double_t t1, t1min, t1max, t2, t2min, t2max, s1, s2, d3;
  Double_t wtreat, ztreat, xintreat[10], xoutreat[10];
  Int_t pdg;

  TLine *line;
  TH1D *h_lpairpp[NHIST], *h_lpairor[NHIST];
  TH1D *htmp;
  TCanvas *c[NHIST];
  Bool_t show[NHIST];
  TLegend *leg;
  TLorentzVector lep1, lep2, prot;
  TPaveText *text;

  for (Int_t i=0; i<NHIST; i++) {
    //if (i<23) show[i] = false;
    //else 
      show[i] = true;
  }
  
  // LPAIR's TTree definition
  Double_t px_[n], py_[n], pz_[n], e_[n], m_[n], eta_[n];
  Int_t pdgId_[n], npart_;
  Int_t insg_;
  Double_t t1_, t1min_, t1max_, t2_, t2min_, t2max_, s1_, s2_, d3_;
  Double_t wtreat_, valtreat_, xtreat_[10], ztreat_[10];
  TFile *lp;
  TTree *tree;
  bool lep1set, lep2set, pset;

  //lp = new TFile("samples/lpair-pt5-mumu-elastic.root");
  //lp = new TFile("samples/lpair-7tev-elastic-pt5.root");
  lp = new TFile("samples/clpair-7tev-elastic-pt5.root");
  //lp = new TFile("samples/lpair-mumu-pt15-8tev-elastic.root");
  //lp = new TFile("samples/lpair-7tev-elastic-nocuts.root");
  //lp = new TFile("samples/lpair-7tev-elastic-pt10-theta0to180.root");
  tree = (TTree*)(lp->Get("h4444"));
  tree->SetBranchAddress("px", px_);
  tree->SetBranchAddress("py", py_);
  tree->SetBranchAddress("pz", pz_);
  tree->SetBranchAddress("E", e_);
  tree->SetBranchAddress("m", m_);
  tree->SetBranchAddress("Eta", eta_);
  tree->SetBranchAddress("icode", pdgId_);
  tree->SetBranchAddress("ip", &npart_);
  tree->SetBranchAddress("t1", &t1_);
  tree->SetBranchAddress("t1min", &t1min_);
  tree->SetBranchAddress("t1max", &t1max_);
  tree->SetBranchAddress("t2", &t2_);
  tree->SetBranchAddress("t2min", &t2min_);
  tree->SetBranchAddress("t2max", &t2max_);
  tree->SetBranchAddress("s1", &s1_);
  tree->SetBranchAddress("s2", &s2_);
  tree->SetBranchAddress("d3", &d3_);
  tree->SetBranchAddress("wtreat", &wtreat_);
  tree->SetBranchAddress("valtreat", &valtreat_);
  tree->SetBranchAddress("xtreat", xtreat_);
  tree->SetBranchAddress("ztreat", ztreat_);
  tree->SetBranchAddress("insetgen", insg_);

  gStyle->SetOptStat(0);

  h_lpairpp[PTSINGLE] = new TH1D("pt", "p_{T}(#mu^{#pm})", 200, 0., 100.);
  h_lpairor[PTSINGLE] = new TH1D("pt_2", "p_{T}(#mu^{#pm})", 200, 0., 100.);
  h_lpairpp[PXSINGLE] = new TH1D("px", "p_{x}(#mu^{#pm})", 200, -100., 100.);
  h_lpairor[PXSINGLE] = new TH1D("px_2", "p_{x}(#mu^{#pm})", 200, -100., 100.);
  h_lpairpp[PYSINGLE] = new TH1D("py", "p_{y}(#mu^{#pm})", 200, -100., 100.);
  h_lpairor[PYSINGLE] = new TH1D("py_2", "p_{y}(#mu^{#pm})", 200, -100., 100.);
  h_lpairpp[PZSINGLE] = new TH1D("pz", "p_{z}(#mu^{#pm})", 200, -100., 100.);
  h_lpairor[PZSINGLE] = new TH1D("pz_2", "p_{z}(#mu^{#pm})", 200, -100., 100.);
  h_lpairpp[ESINGLE] = new TH1D("e", "E (#mu^{#pm})", 200, 0., 100.);
  h_lpairor[ESINGLE] = new TH1D("e_2", "E(#mu^{#pm})", 200, 0., 100.);
  h_lpairpp[PSINGLE] = new TH1D("p", "p(#mu^{#pm})", 200, 0., 100.);
  h_lpairor[PSINGLE] = new TH1D("p_2", "p(#mu^{#pm})", 200, 0., 100.);
  h_lpairpp[ETASINGLE] = new TH1D("eta", "#eta(#mu^{#pm})", 200, -10., 10.);
  h_lpairor[ETASINGLE] = new TH1D("eta_2", "#eta(#mu^{#pm})", 200, -10., 10.);
  h_lpairpp[PHISINGLE] = new TH1D("phi", "#phi(#mu^{#pm})/#pi", 60, -1., 1.);
  h_lpairor[PHISINGLE] = new TH1D("phi_2", "#phi(#mu^{#pm})/#pi", 60, -1., 1.);
  h_lpairpp[THETASINGLE] = new TH1D("theta", "#theta(#mu^{#pm})/#pi", 100, 0., 1.);
  h_lpairor[THETASINGLE] = new TH1D("theta_2", "#theta(#mu^{#pm})/#pi", 100, 0., 1.);
  h_lpairpp[DPT] = new TH1D("dpt", "#Delta p_{T}(#mu^{+}#mu^{-})", 100, 0., 5.);
  h_lpairor[DPT] = new TH1D("dpt_2", "#Delta p_{T}(#mu^{+}#mu^{-})", 100, 0., 5.);
  h_lpairpp[ACOP] = new TH1D("acop", "1-#left|#Delta#phi(#mu^{+}#mu^{-})/#pi#right|", 100, 0., .5);
  h_lpairor[ACOP] = new TH1D("acop_2", "1-#left|#Delta#phi(#mu^{+}#mu^{-})/#pi#right|", 100, 0., .5);
  h_lpairpp[MPAIR] = new TH1D("mass", "m(#mu^{+}#mu^{-})", 200, 0., 100.);
  h_lpairor[MPAIR] = new TH1D("mass_2", "m(#mu^{+}#mu^{-})", 200, 0., 100.);
  h_lpairpp[PTPAIR] = new TH1D("ptpair", "p_{T}(#mu^{+}#mu^{-})", 100, 0., 5.);
  h_lpairor[PTPAIR] = new TH1D("ptpair_2", "p_{T}(#mu^{+}#mu^{-})", 100, 0., 5.);
  h_lpairpp[YPAIR] = new TH1D("ypair", "y(#mu^{+}#mu^{-})", 100, -15., 15.);
  h_lpairor[YPAIR] = new TH1D("ypair_2", "y(#mu^{+}#mu^{-})", 100, -15., 15.);
  h_lpairpp[Q2] = new TH1D("q2m", "Q^{2}", 200, 0., 100.);
  h_lpairor[Q2] = new TH1D("q2m_2", "Q^{2}", 200, 0., 100.);
  h_lpairpp[PPROTON] = new TH1D("pp", "p_{proton}", (Int_t)ep/20, 0., ep);
  h_lpairor[PPROTON] = new TH1D("pp_2", "p_{proton}", (Int_t)ep/20, 0., ep);
  h_lpairpp[T1] = new TH1D("t1", "-t_{1}", 200, 0., 1.);
  h_lpairor[T1] = new TH1D("t1_2", "-t_{1}", 200, 0., 1.);
  h_lpairpp[T1MIN] = new TH1D("t1min", "-t_{1}^{min}", 200, 0., 1.e-2);
  h_lpairor[T1MIN] = new TH1D("t1min_2", "-t_{1}^{min}", 200, 0., 1.e-2);
  h_lpairpp[T1MAX] = new TH1D("t1max", "-t_{1}^{max}", 20, .999e5, 1.001e5);
  h_lpairor[T1MAX] = new TH1D("t1max_2", "-t_{1}^{max}", 20, .999e5, 1.001e5);
  h_lpairpp[T2] = new TH1D("t2", "-t_{2}", 200, 0., 1.);
  h_lpairor[T2] = new TH1D("t2_2", "-t_{2}", 200, 0., 1.);
  h_lpairpp[S1] = new TH1D("s1", "s_{1}", 250, 0., .5e6);
  h_lpairor[S1] = new TH1D("s1_2", "s_{1}", 250, 0., .5e6);
  h_lpairpp[S2] = new TH1D("s2", "s_{2}", 250, 0., .5e6);
  h_lpairor[S2] = new TH1D("s2_2", "s_{2}", 250, 0., .5e6);
  h_lpairpp[D3] = new TH1D("d3", "#delta_{3}", 200, 0., 1.e6);
  h_lpairor[D3] = new TH1D("d3_2", "#delta_{3}", 200, 0., 1.e6);
  h_lpairpp[WTREAT] = new TH1D("wtrt", "w_{treat}", 100, 0., 10);
  h_lpairor[WTREAT] = new TH1D("wtrt_2", "w_{treat}", 100, 0., 10);
  h_lpairpp[ZTREAT] = new TH1D("ztrt", "z_{treat}", 100, 0., 200);
  h_lpairor[ZTREAT] = new TH1D("ztrt_2", "z_{treat}", 100, 0., 200);
  h_lpairpp[XIN0TREAT] = new TH1D("xintrt0", "x^{in}_{treat}[0]", 100, 0., 1.);
  h_lpairor[XIN0TREAT] = new TH1D("xintrt0_2", "x^{in}_{treat}[0]", 100, 0., 1.);
  h_lpairpp[XIN1TREAT] = new TH1D("xintrt1", "x^{in}_{treat}[1]", 100, 0., 1.);
  h_lpairor[XIN1TREAT] = new TH1D("xintrt1_2", "x^{in}_{treat}[1]", 100, 0., 1.);
  h_lpairpp[XIN2TREAT] = new TH1D("xintrt2", "x^{in}_{treat}[2]", 100, 0., 1.);
  h_lpairor[XIN2TREAT] = new TH1D("xintrt2_2", "x^{in}_{treat}[2]", 100, 0., 1.);
  h_lpairpp[XIN3TREAT] = new TH1D("xintrt3", "x^{in}_{treat}[3]", 100, 0., 1.);
  h_lpairor[XIN3TREAT] = new TH1D("xintrt3_2", "x^{in}_{treat}[3]", 100, 0., 1.);
  h_lpairpp[XIN4TREAT] = new TH1D("xintrt4", "x^{in}_{treat}[4]", 100, 0., 1.);
  h_lpairor[XIN4TREAT] = new TH1D("xintrt4_2", "x^{in}_{treat}[4]", 100, 0., 1.);
  h_lpairpp[XIN5TREAT] = new TH1D("xintrt5", "x^{in}_{treat}[5]", 100, 0., 1.);
  h_lpairor[XIN5TREAT] = new TH1D("xintrt5_2", "x^{in}_{treat}[5]", 100, 0., 1.);
  h_lpairpp[XIN6TREAT] = new TH1D("xintrt6", "x^{in}_{treat}[6]", 100, 0., 1.);
  h_lpairor[XIN6TREAT] = new TH1D("xintrt6_2", "x^{in}_{treat}[6]", 100, 0., 1.);
  h_lpairpp[XOUT0TREAT] = new TH1D("xoutrt0", "x^{out}_{treat}[0]", 100, 0., 1.);
  h_lpairor[XOUT0TREAT] = new TH1D("xoutrt0_2", "x^{out}_{treat}[0]", 100, 0., 1.);
  h_lpairpp[XOUT1TREAT] = new TH1D("xoutrt1", "x^{out}_{treat}[1]", 100, 0., 1.);
  h_lpairor[XOUT1TREAT] = new TH1D("xoutrt1_2", "x^{out}_{treat}[1]", 100, 0., 1.);
  h_lpairpp[XOUT2TREAT] = new TH1D("xoutrt2", "x^{out}_{treat}[2]", 100, 0., 1.);
  h_lpairor[XOUT2TREAT] = new TH1D("xoutrt2_2", "x^{out}_{treat}[2]", 100, 0., 1.);
  h_lpairpp[XOUT3TREAT] = new TH1D("xoutrt3", "x^{out}_{treat}[3]", 100, 0., 1.);
  h_lpairor[XOUT3TREAT] = new TH1D("xoutrt3_2", "x^{out}_{treat}[3]", 100, 0., 1.);
  h_lpairpp[XOUT4TREAT] = new TH1D("xoutrt4", "x^{out}_{treat}[4]", 100, 0., 1.);
  h_lpairor[XOUT4TREAT] = new TH1D("xoutrt4_2", "x^{out}_{treat}[4]", 100, 0., 1.);
  h_lpairpp[XOUT5TREAT] = new TH1D("xoutrt5", "x^{out}_{treat}[5]", 100, 0., 1.);
  h_lpairor[XOUT5TREAT] = new TH1D("xoutrt5_2", "x^{out}_{treat}[5]", 100, 0., 1.);
  h_lpairpp[XOUT6TREAT] = new TH1D("xoutrt6", "x^{out}_{treat}[6]", 100, 0., 1.);
  h_lpairor[XOUT6TREAT] = new TH1D("xoutrt6_2", "x^{out}_{treat}[6]", 100, 0., 1.);

  // First fetch the LPAIR++ output
  in.open("test");
  i = 0;
  lep1set = lep2set = false;
  while(in >> e >> px >> py >> pz >> pt >> m >> eta >> pdg >> weight) {
    if (maxEvts>0 && i/2>maxEvts) break;
    if (i%2==0 && (i/2)%10000==0) {
      cout << "[LPAIR++] Event #" << i/2 << endl;
    }
    if (i<5) cout << i << "\t" << pdg << "\t" << m << "\t" << eta << "\t" << px << "\t" << py << "\t" << pz << "\t" << pt << "\t" << e << endl;
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
      h_lpairpp[PTSINGLE]->Fill(lep1.Pt());
      h_lpairpp[PXSINGLE]->Fill(lep1.Px());
      h_lpairpp[PYSINGLE]->Fill(lep1.Py());
      h_lpairpp[PZSINGLE]->Fill(lep1.Pz());
      h_lpairpp[ESINGLE]->Fill(lep1.E());
      h_lpairpp[PSINGLE]->Fill(lep1.P());
      h_lpairpp[ETASINGLE]->Fill(eta);
      h_lpairpp[PHISINGLE]->Fill(lep1.Phi()/pi);
      h_lpairpp[THETASINGLE]->Fill(lep1.Theta()/pi);
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
      h_lpairpp[ACOP]->Fill(1-fabs((lep1.Phi()-lep2.Phi()))/pi);
      h_lpairpp[DPT]->Fill(fabs(lep1.Pt()-lep2.Pt()));
      h_lpairpp[MPAIR]->Fill((lep1+lep2).M());
      h_lpairpp[PTPAIR]->Fill((lep1+lep2).Pt());
      h_lpairpp[YPAIR]->Fill((lep1+lep2).Rapidity());
      lep1set = lep2set = false;
    }
    i++;
  }
  in.close();

  in.open("test_q2");
  while (in >> q2m >> pp3 >> pp5 >> t1 >> t1min >> t1max >> t2 >> t2min >> t2max >> s1 >> s2 >> d3) {
    h_lpairpp[Q2]->Fill(-q2m);
    h_lpairpp[PPROTON]->Fill(pp3);
    h_lpairpp[PPROTON]->Fill(pp5);
    h_lpairpp[T1]->Fill(-t1);
    h_lpairpp[T1MIN]->Fill(-t1min);
    h_lpairpp[T1MAX]->Fill(-t1max);
    h_lpairpp[T2]->Fill(-t2);
    h_lpairpp[S1]->Fill(s1);
    h_lpairpp[S2]->Fill(s2);
    h_lpairpp[D3]->Fill(d3);
  }
  in.close();

  in.open("test_vegas");
  while (in >> wtreat 
	 >> ztreat 
	 >> xoutreat[0] >> xoutreat[1] >> xoutreat[2] >> xoutreat[3] >> xoutreat[4] >> xoutreat[5] >> xoutreat[6]
	 >> xintreat[0] >> xintreat[1] >> xintreat[2] >> xintreat[3] >> xintreat[4] >> xintreat[5] >> xintreat[6]
	 ) {
    h_lpairpp[WTREAT]->Fill(wtreat);
    h_lpairpp[ZTREAT]->Fill(ztreat);
    h_lpairpp[XIN0TREAT]->Fill(xintreat[0]);
    h_lpairpp[XIN1TREAT]->Fill(xintreat[1]);
    h_lpairpp[XIN2TREAT]->Fill(xintreat[2]);
    h_lpairpp[XIN3TREAT]->Fill(xintreat[3]);
    h_lpairpp[XIN4TREAT]->Fill(xintreat[4]);
    h_lpairpp[XIN5TREAT]->Fill(xintreat[5]);
    h_lpairpp[XIN6TREAT]->Fill(xintreat[6]);
    h_lpairpp[XOUT0TREAT]->Fill(xoutreat[0]);
    h_lpairpp[XOUT1TREAT]->Fill(xoutreat[1]);
    h_lpairpp[XOUT2TREAT]->Fill(xoutreat[2]);
    h_lpairpp[XOUT3TREAT]->Fill(xoutreat[3]);
    h_lpairpp[XOUT4TREAT]->Fill(xoutreat[4]);
    h_lpairpp[XOUT5TREAT]->Fill(xoutreat[5]);
    h_lpairpp[XOUT6TREAT]->Fill(xoutreat[6]);
  }
  in.close();
  
  // Then fetch the LPAIR output (converted as a TTree)
  for (i=0; i<tree->GetEntries(); i++) {
    if (maxEvts>0 && i>maxEvts) break;
    if (i%10000==0) {
      cout << "[ LPAIR ] Event #" << i << endl;
    }
    h_lpairor[T1]->Fill(-t1_);
    h_lpairor[T1MIN]->Fill(-t1min_);
    h_lpairor[T1MAX]->Fill(-t1max_);
    h_lpairor[T2]->Fill(-t2_);
    h_lpairor[S1]->Fill(s1_);
    h_lpairor[S2]->Fill(s2_);
    h_lpairor[D3]->Fill(d3_);
    h_lpairor[WTREAT]->Fill(wtreat_);
    h_lpairor[ZTREAT]->Fill(valtreat_);
    h_lpairor[XIN0TREAT]->Fill(xtreat_[0]);
    h_lpairor[XIN1TREAT]->Fill(xtreat_[1]);
    h_lpairor[XIN2TREAT]->Fill(xtreat_[2]);
    h_lpairor[XIN3TREAT]->Fill(xtreat_[3]);
    h_lpairor[XIN4TREAT]->Fill(xtreat_[4]);
    h_lpairor[XIN5TREAT]->Fill(xtreat_[5]);
    h_lpairor[XIN6TREAT]->Fill(xtreat_[6]);
    h_lpairor[XOUT0TREAT]->Fill(ztreat_[0]);
    h_lpairor[XOUT1TREAT]->Fill(ztreat_[1]);
    h_lpairor[XOUT2TREAT]->Fill(ztreat_[2]);
    h_lpairor[XOUT3TREAT]->Fill(ztreat_[3]);
    h_lpairor[XOUT4TREAT]->Fill(ztreat_[4]);
    h_lpairor[XOUT5TREAT]->Fill(ztreat_[5]);
    h_lpairor[XOUT6TREAT]->Fill(ztreat_[6]);
    
    lep1set = lep2set = false;
    pset = false;
    tree->GetEntry(i);
    for (j=0; j<npart_; j++) {
      if (pdgId_[j]==2212) {
        prot.SetXYZM(px_[j], py_[j], pz_[j], m_[j]);
        h_lpairor[PPROTON]->Fill(prot.P());
      }
      if (pdgId_[j]==2212 && !pset) {
        prot.SetXYZM(px_[j], py_[j], pz_[j], m_[j]);
        q2m = -(prot.P()-ep);
        h_lpairor[Q2]->Fill(q2m);
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
      //cout << "filling!" << endl;
      //if (fabs((lep1+lep2).Rapidity())>5.) continue;
      h_lpairor[PTSINGLE]->Fill(lep1.Pt());
      h_lpairor[PXSINGLE]->Fill(lep1.Px());
      h_lpairor[PYSINGLE]->Fill(lep1.Py());
      h_lpairor[PZSINGLE]->Fill(lep1.Pz());
      h_lpairor[ESINGLE]->Fill(lep1.E());
      h_lpairor[PSINGLE]->Fill(lep1.P());
      h_lpairor[ETASINGLE]->Fill(lep1.Eta());
      h_lpairor[PHISINGLE]->Fill(lep1.Phi()/pi);
      h_lpairor[THETASINGLE]->Fill(lep1.Theta()/pi);
      h_lpairor[ACOP]->Fill(1-fabs((lep1.Phi()-lep2.Phi()))/pi);
      h_lpairor[DPT]->Fill(fabs(lep1.Pt()-lep2.Pt()));
      h_lpairor[MPAIR]->Fill((lep1+lep2).M());
      h_lpairor[PTPAIR]->Fill((lep1+lep2).Pt());
      h_lpairor[YPAIR]->Fill((lep1+lep2).Rapidity());
      lep1set = lep2set = false;
    }
  }

  leg = new TLegend(.76, .71, .89, .81);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kBlack);

  text = new TPaveText(.66, .83, .91, .89, "NDC");
  //text = new TPaveText(.65, 1.04, .9, 1.09, "NDC");
  ss.str(""); ss << "LPAIR/LPAIR++ with " << maxEvts << " events";
  text->AddText(ss.str().c_str());
  text->SetFillColor(kWhite);
  text->SetLineColor(kWhite);
  text->SetLineWidth(0);
  text->SetShadowColor(kWhite);
  text->SetTextFont(42);

  Int_t n = 0.1;

  Double_t max;

  for (i=0; i<NHIST; i++) {
    if (!show[i]) continue;
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
    c_2->SetGrid(1,1);

    c[i]->cd(1);
    h_lpairpp[i]->SetFillColor(kBlue);
    //h_lpairpp[i]->SetFillStyle(3005);
    h_lpairpp[i]->SetFillStyle(3003);
    h_lpairpp[i]->SetLineColor(kBlack);
    h_lpairpp[i]->SetLineWidth(1);
    //ss.str(""); ss << "#frac{1}{#sigma} #frac{d#sigma}{d" << h_lpairpp[i]->GetTitle() << "}";
    ss.str(""); ss << "#frac{dN}{d" << h_lpairpp[i]->GetTitle() << "}";
    h_lpairpp[i]->GetXaxis()->SetTitleFont(43);
    h_lpairpp[i]->GetXaxis()->SetTitleSize(14);
    h_lpairpp[i]->GetXaxis()->SetTitleOffset(4.);
    h_lpairpp[i]->GetYaxis()->SetTitleFont(43);
    h_lpairpp[i]->GetYaxis()->SetTitleSize(14);
    h_lpairpp[i]->GetYaxis()->SetTitleOffset(1.4);
    h_lpairpp[i]->GetXaxis()->SetLabelFont(43);
    h_lpairpp[i]->GetXaxis()->SetLabelSize(14);
    h_lpairpp[i]->GetYaxis()->SetLabelFont(43);
    h_lpairpp[i]->GetYaxis()->SetLabelSize(14);
    //h_lpairpp[i]->Scale(1./h_lpairpp[i]->Integral());
    h_lpairor[i]->SetTitle("");
    h_lpairor[i]->GetYaxis()->SetTitle(ss.str().c_str());
    h_lpairor[i]->GetYaxis()->SetTitleFont(43);
    h_lpairor[i]->GetYaxis()->SetTitleSize(14);
    h_lpairor[i]->GetYaxis()->SetTitleOffset(1.4);
    h_lpairor[i]->SetFillColor(kRed);
    //h_lpairor[i]->SetFillStyle(3004);
    h_lpairor[i]->SetFillStyle(3001);
    h_lpairor[i]->SetLineColor(kBlack);
    h_lpairor[i]->SetLineWidth(1);

    //h_lpairor[i]->Scale(1./h_lpairor[i]->Integral());
    h_lpairor[i]->Draw();
    h_lpairpp[i]->Draw("SAME");
    max = TMath::Max(h_lpairor[i]->GetBinContent(h_lpairor[i]->GetMaximumBin()), h_lpairpp[i]->GetBinContent(h_lpairpp[i]->GetMaximumBin()));
    h_lpairor[i]->GetYaxis()->SetRangeUser(.01, max*1.2);
    if (n==0) {
      leg->AddEntry(h_lpairpp[i], "LPAIR++");
      leg->AddEntry(h_lpairor[i], "LPAIR");
    }
    leg->Draw("SAME");
    text->Draw();

    c[i]->cd(2);
    line = new TLine(h_lpairpp[i]->GetXaxis()->GetXmin(),1., h_lpairpp[i]->GetXaxis()->GetXmax(), 1.);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    htmp = (TH1D*)h_lpairpp[i]->Clone();
    htmp->Divide(h_lpairor[i]);
    htmp->SetTitle("");
    htmp->GetXaxis()->SetTitle(h_lpairpp[i]->GetTitle());
    htmp->GetYaxis()->SetTitle("LPAIR++/LPAIR");
    htmp->Draw("E");
    line->Draw();
    n++;
  }
}
