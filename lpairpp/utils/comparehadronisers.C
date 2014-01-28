#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLorentzVector.h"

#define NREMN 0
#define REMN_MX 1
#define REMN_ETA 2
#define REMN_PHI 3
#define REMN_PX 4
#define REMN_PY 5
#define REMN_PZ 6
#define REMN_PT 7
#define NHIST 8

void comparehadronisers() {

  const Int_t n = 5000;

  TFile *f_js, *f_py;
  TTree *t_js, *t_py;

  Int_t np_js;
  Double_t eta_js[n], phi_js[n], rapidity_js[n], charge_js[n];
  Double_t px_js[n], py_js[n], pz_js[n], pt_js[n], E_js[n], M_js[n];
  Int_t PID_js[n], role_js[n], parentid_js[n], isstable_js[n];

  Int_t np_py;
  Double_t eta_py[n], phi_py[n], rapidity_py[n], charge_py[n];
  Double_t px_py[n], py_py[n], pz_py[n], pt_py[n], E_py[n], M_py[n];
  Int_t PID_py[n], role_py[n], parentid_py[n], isstable_py[n];

  TH1D *h_js[NHIST], *h_py[NHIST];

  TCanvas *c[NHIST];
  TLegend *l[NHIST];

  TLorentzVector part;

  TPaveText *text;
  std::stringstream ss;

  Int_t num_remn;

  f_js = new TFile("events_lpairpp_jetset_100kevts.root");
  f_py = new TFile("events_lpairpp_pythia_100kevts.root");

  t_js = (TTree*)f_js->Get("h4444");
  t_py = (TTree*)f_py->Get("h4444");

  TString title[NHIST];
  title[NREMN] = "Number of particles in the proton remnants";
  title[REMN_MX] = "M_{X}";
  title[REMN_ETA] = "#eta^{remnants}";
  title[REMN_PHI] = "#phi^{remnants}";
  title[REMN_PX] = "p_{x}^{remnants}";
  title[REMN_PY] = "p_{y}^{remnants}";
  title[REMN_PZ] = "p_{z}^{remnants}";
  title[REMN_PT] = "p_{T}^{remnants}";
    
  h_js[NREMN] = new TH1D("h_remn_js", "", 50, -.5, 49.5);
  h_py[NREMN] = new TH1D("h_remn_py", "", 50, -.5, 49.5);
  h_js[REMN_MX] = new TH1D("h_mtot_remn_js", "", 200, 0., 50);
  h_py[REMN_MX] = new TH1D("h_mtot_remn_py", "", 200, 0., 50);
  h_js[REMN_ETA] = new TH1D("h_eta_remn_js", "", 120, -15., 15.);
  h_py[REMN_ETA] = new TH1D("h_eta_remn_py", "", 120, -15., 15.);
  h_js[REMN_PHI] = new TH1D("h_phi_remn_js", "", 20, -5., 5.);
  h_py[REMN_PHI] = new TH1D("h_phi_remn_py", "", 20, -5., 5.);
  h_js[REMN_PX] = new TH1D("h_px_remn_js", "", 200, -5., 5.);
  h_py[REMN_PX] = new TH1D("h_px_remn_py", "", 200, -5., 5.);
  h_js[REMN_PY] = new TH1D("h_py_remn_js", "", 200, -5., 5.);
  h_py[REMN_PY] = new TH1D("h_py_remn_py", "", 200, -5., 5.);
  h_js[REMN_PZ] = new TH1D("h_pz_remn_js", "", 350, 0., 3500.);
  h_py[REMN_PZ] = new TH1D("h_pz_remn_py", "", 350, 0., 3500.);
  h_js[REMN_PT] = new TH1D("h_pt_remn_js", "", 200, 0., 50.);
  h_py[REMN_PT] = new TH1D("h_pt_remn_py", "", 200, 0., 50.);

  t_js->SetBranchAddress("npart", &np_js);
  t_js->SetBranchAddress("Eta", eta_js);
  t_js->SetBranchAddress("phi", phi_js);
  t_js->SetBranchAddress("rapidity", rapidity_js);
  t_js->SetBranchAddress("px", px_js);
  t_js->SetBranchAddress("py", py_js);
  t_js->SetBranchAddress("pz", pz_js);
  t_js->SetBranchAddress("pt", pt_js);
  t_js->SetBranchAddress("icode", PID_js);
  t_js->SetBranchAddress("role", role_js);
  t_js->SetBranchAddress("parent", parentid_js);
  t_js->SetBranchAddress("stable", isstable_js);
  t_js->SetBranchAddress("E", E_js);
  t_js->SetBranchAddress("m", M_js);
  t_js->SetBranchAddress("charge", charge_js);

  t_py->SetBranchAddress("npart", &np_py);
  t_py->SetBranchAddress("Eta", eta_py);
  t_py->SetBranchAddress("phi", phi_py);
  t_py->SetBranchAddress("rapidity", rapidity_py);
  t_py->SetBranchAddress("px", px_py);
  t_py->SetBranchAddress("py", py_py);
  t_py->SetBranchAddress("pz", pz_py);
  t_py->SetBranchAddress("pt", pt_py);
  t_py->SetBranchAddress("icode", PID_py);
  t_py->SetBranchAddress("role", role_py);
  t_py->SetBranchAddress("parent", parentid_py);
  t_py->SetBranchAddress("stable", isstable_py);
  t_py->SetBranchAddress("E", E_py);
  t_py->SetBranchAddress("m", M_py);
  t_py->SetBranchAddress("charge", charge_py);

  for (int e=0; e<t_js->GetEntries(); e++) {
    if (e%10000==0) std::cout << "--> " << e << std::endl;
    t_js->GetEntry(e);
    num_remn = 0;
    TLorentzVector remn;
    for (int p=0; p<np_js; p++) {
      //std::cout << p << "-> " << (role_js[p]==3&&isstable_js[p]) << std::endl;
      if (role_js[p]==3 && isstable_js[p]) {
	h_js[REMN_ETA]->Fill(eta_js[p]);
	h_js[REMN_PHI]->Fill(phi_js[p]);
	h_js[REMN_PT]->Fill(pt_js[p]);
	h_js[REMN_PX]->Fill(px_js[p]);
	h_js[REMN_PY]->Fill(py_js[p]);
	h_js[REMN_PZ]->Fill(pz_js[p]);
	part.SetXYZM(px_js[p], py_js[p], pz_js[p], M_js[p]);
	remn += part;
	num_remn++;
      }
    }
    h_js[REMN_MX]->Fill(remn.M());
    h_js[NREMN]->Fill(num_remn-.5);
  }

  for (int e=0; e<t_py->GetEntries(); e++) {
    if (e%10000==0) std::cout << "--> " << e << std::endl;
    t_py->GetEntry(e);
    num_remn = 0;
    TLorentzVector remn;
    for (int p=0; p<np_py; p++) {
      //std::cout << p << "-> " << (role_py[p]==3&&isstable_py[p]) << std::endl;
      if (role_py[p]==3 && isstable_py[p]) {
	h_py[REMN_ETA]->Fill(eta_py[p]);
	h_py[REMN_PHI]->Fill(phi_py[p]);
	h_py[REMN_PT]->Fill(pt_py[p]);
	h_py[REMN_PX]->Fill(px_py[p]);
	h_py[REMN_PY]->Fill(py_py[p]);
	h_py[REMN_PZ]->Fill(pz_py[p]);
	part.SetXYZM(px_py[p], py_py[p], pz_py[p], M_py[p]);
	remn += part;
	num_remn++;
      }
    }
    h_py[REMN_MX]->Fill(remn.M());
    h_py[NREMN]->Fill(num_remn-.5);
  }

  gStyle->SetOptStat(0);

  for (int p=0; p<NHIST; p++) {

    text = new TPaveText(.4, .92, .93, .96, "NDC");
    text->SetTextAlign(33);
    ss.str(""); ss << "LPAIR++ with " << t_js->GetEntries() << " events";
    text->AddText(ss.str().c_str());
    text->SetFillColor(kWhite);
    text->SetLineColor(kWhite);
    text->SetLineWidth(0);
    text->SetShadowColor(kWhite);
    text->SetTextFont(42);

    c[p] = new TCanvas(h_js[p]->GetName());
    l[p] = new TLegend(.6, .72, .78, .82);
    l[p]->SetFillColor(kWhite);
    l[p]->SetLineColor(kWhite);
    l[p]->SetTextFont(42);
    h_js[p]->SetLineColor(kRed);
    h_js[p]->Draw();
    h_js[p]->GetXaxis()->SetTitle(title[p]);
    h_py[p]->Draw("SAME");
    l[p]->AddEntry(h_js[p], "Jetset 7.410");
    l[p]->AddEntry(h_py[p], "Pythia 6.4.28");
    l[p]->Draw("SAME");
    text->Draw();
  }

}
