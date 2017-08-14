#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TLorentzVector.h"

#include "TStyle.h"

#define PT_SINGLEL 0
#define PX_SINGLEL 1
#define PY_SINGLEL 2
#define PZ_SINGLEL 3
#define E_SINGLEL 4
#define P_SINGLEL 5
#define ETA_SINGLEL 6
#define PHI_SINGLEL 7
#define THETA_SINGLEL 8
#define DPT_DIL 9
#define ACOP_DIL 10
#define M_DIL 11
#define PT_DIL 12
#define RAP_DIL 13
#define P_PROTON 14
#define REMN_NUM 15
#define REMN_NUM_CH 16
#define REMN_NUM_NT 17
#define REMN_TOT_M 18
#define REMN_TOT_M_HAD 19
#define REMN_P 20
#define REMN_PT 21
#define REMN_E 22

#define NHIST 23

using namespace std;

void distributions()
{
  //SetStyle();

  const Int_t lepPdg = 13;
  const Int_t n=500;
  //const Int_t maxEvts = -1;
  const Int_t maxEvts = 1e5;
  const Double_t ep = 3500.;
  const Double_t pi = acos(-1);

  Double_t max;
  TPad *c_1, *c_2;

  stringstream ss;

  TLine *line;
  TH1D *h_lpairpp[NHIST], *h_lpairor[NHIST];
  TH1D *htmp;
  TCanvas *c[NHIST];
  Bool_t show[NHIST];
  TLegend *leg;
  TString hadroniser;
  TLorentzVector lep1, lep2, remn, tot_remn, prot;
  TPaveText *text;

  for (i=0; i<NHIST; i++) {
    //if (i==24 || (i>=25 && i<=31)) show[i] = false;
    //else 
      show[i] = true;
  }
  Int_t nremn, nremn_ch, nremn_nt;
  
  // LPAIR's TTree definition
  Double_t xsect_lp, errxsect_lp, mx_lp;
  Double_t txsect_lp, terrxsect_lp;
  Double_t px_lp[n], py_lp[n], pz_lp[n], e_lp[n], m_lp[n], eta_lp[n];
  Int_t pdgId_lp[n], npart_lp, stable_lp[n], mother_lp[n];
  Double_t charge_lp[n];

  Int_t n_lp;
  TFile *lp, *clp;
  TTree *t_lp, *t_clp;
  bool lep1set, lep2set, pset;

  //lp = new TFile("samples/clpair-7tev-elastic-pt5.root");
  lp = new TFile("samples/lpair-7tev-singlediss-pt5.root");
  t_lp = (TTree*)(lp->Get("h4444"));
  t_lp->SetBranchAddress("ip", &npart_lp);
  t_lp->SetBranchAddress("xsect", &txsect_lp);
  t_lp->SetBranchAddress("errxsect", &terrxsect_lp);
  t_lp->SetBranchAddress("MX", &mx_lp);
  t_lp->SetBranchAddress("px", px_lp);
  t_lp->SetBranchAddress("py", py_lp);
  t_lp->SetBranchAddress("pz", pz_lp);
  t_lp->SetBranchAddress("E", e_lp);
  t_lp->SetBranchAddress("m", m_lp);
  t_lp->SetBranchAddress("Eta", eta_lp);
  t_lp->SetBranchAddress("stable", stable_lp);
  t_lp->SetBranchAddress("icode", pdgId_lp);
  t_lp->SetBranchAddress("charge", charge_lp);
  t_lp->SetBranchAddress("parent", mother_lp);

  hadroniser = "Jetset 7.410";
  //hadroniser = "Pythia 6.4.28";

  //clp = new TFile("events_lpairpp_jetset_100kevts.root");
  //clp = new TFile("events_lpairpp_pythia_100kevts.root");
  clp = new TFile("events_lpairpp_jetset.root");
  //clp = new TFile("events_lpairpp_pythia.root");
  t_clp = (TTree*)(clp->Get("h4444"));
  CepGen::TreeEvent ev;
  ev.attach( t_lp );

  gStyle->SetOptStat(0);

  h_lpairpp[PT_SINGLEL] = new TH1D("pt", "p_{T}(#mu^{#pm})", 200, 0., 100.);
  h_lpairor[PT_SINGLEL] = new TH1D("pt_2", "p_{T}(#mu^{#pm})", 200, 0., 100.);
  h_lpairpp[PX_SINGLEL] = new TH1D("px", "p_{x}(#mu^{#pm})", 200, -100., 100.);
  h_lpairor[PX_SINGLEL] = new TH1D("px_2", "p_{x}(#mu^{#pm})", 200, -100., 100.);
  h_lpairpp[PY_SINGLEL] = new TH1D("py", "p_{y}(#mu^{#pm})", 200, -100., 100.);
  h_lpairor[PY_SINGLEL] = new TH1D("py_2", "p_{y}(#mu^{#pm})", 200, -100., 100.);
  h_lpairpp[PZ_SINGLEL] = new TH1D("pz", "p_{z}(#mu^{#pm})", 200, -100., 100.);
  h_lpairor[PZ_SINGLEL] = new TH1D("pz_2", "p_{z}(#mu^{#pm})", 200, -100., 100.);
  h_lpairpp[E_SINGLEL] = new TH1D("e", "E (#mu^{#pm})", 200, 0., 100.);
  h_lpairor[E_SINGLEL] = new TH1D("e_2", "E(#mu^{#pm})", 200, 0., 100.);
  h_lpairpp[P_SINGLEL] = new TH1D("p", "p(#mu^{#pm})", 200, 0., 100.);
  h_lpairor[P_SINGLEL] = new TH1D("p_2", "p(#mu^{#pm})", 200, 0., 100.);
  h_lpairpp[ETA_SINGLEL] = new TH1D("eta", "#eta(#mu^{#pm})", 200, -10., 10.);
  h_lpairor[ETA_SINGLEL] = new TH1D("eta_2", "#eta(#mu^{#pm})", 200, -10., 10.);
  h_lpairpp[PHI_SINGLEL] = new TH1D("phi", "#phi(#mu^{#pm})/#pi", 60, -1., 1.);
  h_lpairor[PHI_SINGLEL] = new TH1D("phi_2", "#phi(#mu^{#pm})/#pi", 60, -1., 1.);
  h_lpairpp[THETA_SINGLEL] = new TH1D("theta", "#theta(#mu^{#pm})/#pi", 100, 0., 1.);
  h_lpairor[THETA_SINGLEL] = new TH1D("theta_2", "#theta(#mu^{#pm})/#pi", 100, 0., 1.);
  h_lpairpp[DPT_DIL] = new TH1D("dpt", "#Delta p_{T}(#mu^{+}#mu^{-})", 100, 0., 5.);
  h_lpairor[DPT_DIL] = new TH1D("dpt_2", "#Delta p_{T}(#mu^{+}#mu^{-})", 100, 0., 5.);
  h_lpairpp[ACOP_DIL] = new TH1D("acop", "1-#left|#Delta#phi(#mu^{+}#mu^{-})/#pi#right|", 100, 0., .5);
  h_lpairor[ACOP_DIL] = new TH1D("acop_2", "1-#left|#Delta#phi(#mu^{+}#mu^{-})/#pi#right|", 100, 0., .5);
  h_lpairpp[M_DIL] = new TH1D("mass", "m(#mu^{+}#mu^{-})", 200, 0., 100.);
  h_lpairor[M_DIL] = new TH1D("mass_2", "m(#mu^{+}#mu^{-})", 200, 0., 100.);
  h_lpairpp[PT_DIL] = new TH1D("ptpair", "p_{T}(#mu^{+}#mu^{-})", 100, 0., 5.);
  h_lpairor[PT_DIL] = new TH1D("ptpair_2", "p_{T}(#mu^{+}#mu^{-})", 100, 0., 5.);
  h_lpairpp[RAP_DIL] = new TH1D("ypair", "y(#mu^{+}#mu^{-})", 100, -15., 15.);
  h_lpairor[RAP_DIL] = new TH1D("ypair_2", "y(#mu^{+}#mu^{-})", 100, -15., 15.);
  h_lpairpp[P_PROTON] = new TH1D("pp", "p_{proton}", (Int_t)ep/20, 0., ep);
  h_lpairor[P_PROTON] = new TH1D("pp_2", "p_{proton}", (Int_t)ep/20, 0., ep);
  h_lpairpp[REMN_NUM] = new TH1D("rm_num", "N^{remnants}", 60, -.5, 59.5);
  h_lpairor[REMN_NUM] = new TH1D("rm_num_2", "N^{remnants}", 60, -.5, 59.5);
  h_lpairpp[REMN_NUM_CH] = new TH1D("rm_num_ch", "N_{charged}^{remnants}", 60, -.5, 59.5);
  h_lpairor[REMN_NUM_CH] = new TH1D("rm_num_ch_2", "N_{charged}^{remnants}", 60, -.5, 59.5);
  h_lpairpp[REMN_NUM_NT] = new TH1D("rm_num_nt", "N_{neutral}^{remnants}", 60, -.5, 59.5);
  h_lpairor[REMN_NUM_NT] = new TH1D("rm_num_nt_2", "N_{neutral}^{remnants}", 60, -.5, 59.5);
  h_lpairpp[REMN_TOT_M] = new TH1D("rm_tot_mass", "M_{X}", 175, 0., 350.);
  h_lpairor[REMN_TOT_M] = new TH1D("rm_tot_mass_2", "M_{X}", 175, 0., 350.);
  h_lpairpp[REMN_TOT_M_HAD] = new TH1D("rm_tot_mass_had", "M_{X} (hadronised)", 175, 0., 350.);
  h_lpairor[REMN_TOT_M_HAD] = new TH1D("rm_tot_mass_had_2", "M_{X} (hadronised)", 175, 0., 350.);
  h_lpairpp[REMN_P] = new TH1D("rm_p", "p^{remnants}", (Int_t)ep/20, 0., ep);
  h_lpairor[REMN_P] = new TH1D("rm_p_2", "p^{remnants}", (Int_t)ep/20, 0., ep);
  h_lpairpp[REMN_PT] = new TH1D("rm_pt", "p_{T}^{remnants}", 50, 0., 200.);
  h_lpairor[REMN_PT] = new TH1D("rm_pt_2", "p_{T}^{remnants}", 50, 0., 200.);
  h_lpairpp[REMN_E] = new TH1D("rm_e", "E^{remnants}", 175, 0., 3500.);
  h_lpairor[REMN_E] = new TH1D("rm_e_2", "E^{remnants}", 175, 0., 3500.);

  const unsigned long long n_clp = t_clp->GetEntries(), n_lp = t_lp->GetEntries();

  // First fetch the LPAIR++ output
  for ( unsigned long long i=0; i<n_clp; i++) {
    if ( maxEvts > 0 && i > maxEvts ) break;

    tot_remn.SetXYZM(0., 0., 0., 0.);
    lep1set = lep2set = false;
    nremn = nremn_ch = nremn_nt = 0;

    t_clp->GetEntry( i );

    if ( i == 0 ) {
      xsect_clp = ev.xsect;
      errxsect_clp = ev.errxsect;
      cout << "[LPAIR++] Sigma = " << xsect_clp << " +/- " << errxsect_clp << endl;
    }
    else if (i%20000==0) {
      cout << "[LPAIR++] Event #" << i << endl;
    }

    for ( unsigned short j=0; j<ev.np; j++ ) {
      if ( ev.PID[j] == -lepPdg ) {
        lep1.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] );
        lep1set = true;
      }
      else if ( ev.PID[j] == lepPdg ) {
        lep2.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] );
        lep2set = true;
      }
      //else if (( ev.stable[j] == 1 ) && ( ev.mother[j]!=-1) && (mother_clp[j]!=0) && (mother_clp[j]!=1)) {
      else if ( ev.role[j] == 3 && ev.stable[j] ) {
        // Stable proton remnants
        remn.SetPtEtaPhiM( ev.px[j], ev.eta[j], ev.phi[j], ev.M[j] );
        h_lpairpp[REMN_P]->Fill( remn.P() );
        h_lpairpp[REMN_PT]->Fill( remn.Pt() );
        h_lpairpp[REMN_E]->Fill( remn.E() );
        tot_remn += remn;
        if ( ev.charge[j] == 0. ) nremn_nt++;
        else nremn_ch++;
        nremn++;
      }
    }
    if (lep1set && lep2set) {
      h_lpairpp[PT_SINGLEL]->Fill(lep1.Pt());
      h_lpairpp[PX_SINGLEL]->Fill(lep1.Px());
      h_lpairpp[PY_SINGLEL]->Fill(lep1.Py());
      h_lpairpp[PZ_SINGLEL]->Fill(lep1.Pz());
      h_lpairpp[E_SINGLEL]->Fill(lep1.E());
      h_lpairpp[P_SINGLEL]->Fill(lep1.P());
      h_lpairpp[ETA_SINGLEL]->Fill(lep1.Eta());
      h_lpairpp[PHI_SINGLEL]->Fill(lep1.Phi()/pi);
      h_lpairpp[THETA_SINGLEL]->Fill(lep1.Theta()/pi);
      h_lpairpp[ACOP_DIL]->Fill(1-fabs((lep1.Phi()-lep2.Phi()))/pi);
      h_lpairpp[DPT_DIL]->Fill(fabs(lep1.Pt()-lep2.Pt()));
      h_lpairpp[M_DIL]->Fill((lep1+lep2).M());
      h_lpairpp[PT_DIL]->Fill((lep1+lep2).Pt());
      h_lpairpp[RAP_DIL]->Fill((lep1+lep2).Rapidity());
      h_lpairpp[REMN_NUM]->Fill(nremn-.5);
      h_lpairpp[REMN_NUM_CH]->Fill(nremn_ch-.5);
      h_lpairpp[REMN_NUM_NT]->Fill(nremn_nt-.5);
      lep1set = lep2set = false;
    }
    h_lpairpp[REMN_TOT_M]->Fill(mx_clp);
    h_lpairpp[REMN_TOT_M_HAD]->Fill(tot_remn.M());
  }

  // Then fetch the LPAIR output (as a TTree)
  for (i=0; i<n_lp; i++) {
    if (maxEvts>0 && i>maxEvts) break;

    t_lp->GetEntry(i);

    if (i==0) {
      xsect_lp = txsect_lp;
      errxsect_lp = terrxsect_lp;
      cout << "[ LPAIR ] Sigma = " << xsect_lp << " +/- " << errxsect_lp << endl;
    }
    else if (i%20000==0) {
      cout << "[ LPAIR ] Event #" << i << endl;
    }

    tot_remn.SetXYZM(0., 0., 0., 0.);
    lep1set = lep2set = false;
    pset = false;
    nremn = nremn_ch = nremn_nt = 0;

    for (j=0; j<npart_lp; j++) {
      if (pdgId_lp[j]==-lepPdg) {
	lep1.SetXYZM(px_lp[j], py_lp[j], pz_lp[j], m_lp[j]);
	lep1set = true;
      }
      else if (pdgId_lp[j]==lepPdg) {
	lep2.SetXYZM(px_lp[j], py_lp[j], pz_lp[j], m_lp[j]);
	lep2set = true;
      }
      else if ((stable_lp[j]==1) && (mother_lp[j]!=0) && (mother_lp[j]!=2)) {
	// Stable proton remnants
	remn.SetXYZM(px_lp[j], py_lp[j], pz_lp[j], m_lp[j]);
	h_lpairor[REMN_P]->Fill(remn.P());
	h_lpairor[REMN_PT]->Fill(remn.Pt());
	h_lpairor[REMN_E]->Fill(remn.E());
	tot_remn += remn;
	if (charge_lp[j]==0.) nremn_nt++;
	else nremn_ch++;
	nremn++;
      }
    }
    if (lep1set && lep2set) {
      h_lpairor[PT_SINGLEL]->Fill(lep1.Pt());
      h_lpairor[PX_SINGLEL]->Fill(lep1.Px());
      h_lpairor[PY_SINGLEL]->Fill(lep1.Py());
      h_lpairor[PZ_SINGLEL]->Fill(lep1.Pz());
      h_lpairor[E_SINGLEL]->Fill(lep1.E());
      h_lpairor[P_SINGLEL]->Fill(lep1.P());
      h_lpairor[ETA_SINGLEL]->Fill(lep1.Eta());
      h_lpairor[PHI_SINGLEL]->Fill(lep1.Phi()/pi);
      h_lpairor[THETA_SINGLEL]->Fill(lep1.Theta()/pi);
      h_lpairor[ACOP_DIL]->Fill(1-fabs((lep1.Phi()-lep2.Phi()))/pi);
      h_lpairor[DPT_DIL]->Fill(fabs(lep1.Pt()-lep2.Pt()));
      h_lpairor[M_DIL]->Fill((lep1+lep2).M());
      h_lpairor[PT_DIL]->Fill((lep1+lep2).Pt());
      h_lpairor[RAP_DIL]->Fill((lep1+lep2).Rapidity());
      h_lpairor[REMN_NUM]->Fill(nremn-.5);
      h_lpairor[REMN_NUM_CH]->Fill(nremn_ch-.5);
      h_lpairor[REMN_NUM_NT]->Fill(nremn_nt-.5);
      lep1set = lep2set = false;
    }
    h_lpairor[REMN_TOT_M]->Fill(mx_lp);
    h_lpairor[REMN_TOT_M_HAD]->Fill(tot_remn.M());
  }

  leg = new TLegend(.63, .71, .94, .84);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kBlack);

  //text = new TPaveText(.66, .83, .91, .89, "NDC");
  //text = new TPaveText(.65, 1.04, .9, 1.09, "NDC");
  text = new TPaveText(.4, .925, 1., .98, "NDC");
  text->SetTextAlign(33);
  ss.str(""); ss << "LPAIR/LPAIR++ with " << maxEvts << " events";
  text->AddText(ss.str().c_str());
  text->SetFillColor(kWhite);
  text->SetLineColor(kWhite);
  text->SetLineWidth(0);
  text->SetShadowColor(kWhite);
  text->SetTextFont(42);

  for (i=0; i<NHIST; i++) {
    if (!show[i]) continue;
    c[i] = new TCanvas(); 
    c[i]->Divide(1,2);

    c_1 = (TPad*)(c[i]->GetPad(1));
    c_1->SetPad(0.,.250,1.,1.);
    c_1->SetRightMargin(0.03);
    c_1->SetBottomMargin(0.);
    c_1->SetGrid(1,1);
    c_2 = (TPad*)(c[i]->GetPad(2));
    c_2->SetPad(0.,0.,1.,.250);
    c_2->SetBottomMargin(0.3);
    c_2->SetRightMargin(0.03);
    c_2->SetTopMargin(0.);
    c_2->SetGrid(1,1);

    c[i]->cd(1);
    h_lpairpp[i]->Sumw2();
    ss.str(""); ss << "#frac{1}{#sigma} #frac{d#sigma}{d" << h_lpairpp[i]->GetTitle() << "}";
    //ss.str(""); ss << "#frac{dN}{d" << h_lpairpp[i]->GetTitle() << "}";
    h_lpairpp[i]->GetXaxis()->SetTitleFont(43);
    h_lpairpp[i]->GetXaxis()->SetTitleSize(16);
    h_lpairpp[i]->GetXaxis()->SetTitleOffset(4.);
    h_lpairpp[i]->GetYaxis()->SetTitleFont(43);
    h_lpairpp[i]->GetYaxis()->SetTitleSize(16);
    h_lpairpp[i]->GetYaxis()->SetTitleOffset(1.2);
    h_lpairpp[i]->GetXaxis()->SetLabelFont(43);
    h_lpairpp[i]->GetXaxis()->SetLabelSize(18);
    h_lpairpp[i]->GetYaxis()->SetLabelFont(43);
    h_lpairpp[i]->GetYaxis()->SetLabelSize(18);
    h_lpairpp[i]->SetFillColor(kRed);
    //h_lpairpp[i]->SetFillStyle(3005);
    //h_lpairpp[i]->SetFillStyle(3002);
    h_lpairpp[i]->SetFillStyle(3001);
    h_lpairpp[i]->SetLineColor(kBlack);
    h_lpairpp[i]->SetLineWidth(1);
    //h_lpairpp[i]->SetLineColor(kGray+3);
    h_lpairpp[i]->SetLineStyle(1);
    //h_lpairpp[i]->Scale(1./h_lpairpp[i]->Integral());
    h_lpairor[i]->Sumw2();
    h_lpairor[i]->SetTitle("");
    h_lpairor[i]->GetYaxis()->SetTitle(ss.str().c_str());
    h_lpairor[i]->GetYaxis()->SetTitleFont(43);
    h_lpairor[i]->GetYaxis()->SetTitleSize(16);
    h_lpairor[i]->GetYaxis()->SetTitleOffset(1.2);
    h_lpairor[i]->GetXaxis()->SetLabelFont(43);
    h_lpairor[i]->GetXaxis()->SetLabelSize(18);
    h_lpairor[i]->GetYaxis()->SetLabelFont(43);
    h_lpairor[i]->GetYaxis()->SetLabelSize(18);
    h_lpairor[i]->SetFillColor(kBlue-10);
    //h_lpairor[i]->SetFillStyle(3004);
    //h_lpairor[i]->SetFillStyle(3002);
    //h_lpairor[i]->SetFillStyle(3001);
    h_lpairor[i]->SetLineColor(kBlack);
    h_lpairor[i]->SetLineWidth(1);
    h_lpairor[i]->SetLineStyle(2);

    //h_lpairor[i]->Scale(1./h_lpairor[i]->Integral());
    h_lpairor[i]->Draw("HIST");
    h_lpairpp[i]->Draw("HIST SAME");
    if (i==0) {
      TString leg_lpairpp = "LPAIR++";
      if (hadroniser!="") {
	leg_lpairpp+= " ("+hadroniser+")";
      }
      leg->AddEntry(h_lpairor[i], "LPAIR", "F");
      leg->AddEntry(h_lpairpp[i], leg_lpairpp, "F");
    }

    c[i]->cd(2);
    line = new TLine(h_lpairpp[i]->GetXaxis()->GetXmin(),1., h_lpairpp[i]->GetXaxis()->GetXmax(), 1.);
    //line->SetLineColor(kBlue);
    //line->SetLineWidth(1);
    htmp = (TH1D*)h_lpairpp[i]->Clone();
    htmp->SetFillStyle(3001);
    htmp->SetFillColor(kBlue-10);
    htmp->SetMarkerStyle(8);
    htmp->SetMarkerSize(.6);
    htmp->Divide((TH1D*)h_lpairpp[i]->Clone(), (TH1D*)h_lpairor[i]->Clone(), 1, 1, "B");
    htmp->SetTitle("");
    htmp->GetXaxis()->SetTitle(h_lpairpp[i]->GetTitle());
    htmp->GetYaxis()->SetTitle("LPAIR++/LPAIR");
    htmp->GetYaxis()->SetLabelSize(14);
    htmp->GetYaxis()->SetTitleOffset(1.4);
    htmp->Draw("E");
    line->Draw();
    htmp->Draw("E3 SAME");
    htmp->Draw("E SAME");

    c[i]->cd(1);
    h_lpairpp[i]->Scale(xsect_clp/n_clp);
    h_lpairor[i]->Scale(xsect_lp/n_lp);
    //max = TMath::Max(h_lpairor[i]->GetBinContent(h_lpairor[i]->GetMaximumBin()), h_lpairpp[i]->GetBinContent(h_lpairpp[i]->GetMaximumBin()));
    //h_lpairor[i]->GetYaxis()->SetRangeUser(.01, max*1.2);
    leg->Draw("SAME");
    text->Draw();

  }
}
