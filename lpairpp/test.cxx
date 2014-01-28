#include <iostream>

#include "include/mcgen.h"

// ROOT includes
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

/**
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * @date 27 jan 2014
 */
int main() {
  Parameters ip;
  Event ev;
  double xsec, err;
  Particles particles;
  Particles::iterator p;
  //TH1D *h_npart_outgoing, *h_npart_outgoing_charged, *h_npart_outgoing_neutral;
  TCanvas *c_npart;
  TLegend *leg;
  TTree *tree;
  //Pythia6Hadroniser had;
  Jetset7Hadroniser had;

  ip.in1p = 3500.;
  ip.in2p = 3500.;
  ip.pair = 13;
  ip.p1mod = 11;
  ip.p2mod = 2;
  ip.mcut = 2;
  ip.minenergy = 0.; //FIXME
  ip.minpt = 5.;
  ip.maxgen = 1e1;
  ip.hadroniser = &had;
  //ip.ncvg = 5e3; //FIXME
  //ip.maxgen = 1e5;
  //ip.SetEtaRange(-2.5, 2.5);

  ip.generation = true;
  ip.Dump();

  MCGen mg(&ip);

  mg.ComputeXsection(&xsec, &err);

  const int maxpart = 1000;

  int np;
  double eta[maxpart], phi[maxpart], rapidity[maxpart];
  double px[maxpart], py[maxpart], pz[maxpart], pt[maxpart], E[maxpart], M[maxpart], charge[maxpart];
  int PID[maxpart], parentid[maxpart], isstable[maxpart], role[maxpart];

  tree = new TTree("h4444", "A TTree containing information from the events produced from LPAIR++");
  tree->Branch("npart", &np, "npart/I");
  tree->Branch("Eta", eta, "eta[npart]/D");
  tree->Branch("phi", phi, "phi[npart]/D");
  tree->Branch("rapidity", rapidity, "rapidity[npart]/D");
  tree->Branch("px", px, "px[npart]/D");
  tree->Branch("py", py, "py[npart]/D");
  tree->Branch("pz", pz, "pz[npart]/D");
  tree->Branch("pt", pt, "pt[npart]/D");
  tree->Branch("icode", PID, "PID[npart]/I");
  tree->Branch("role", role, "role[npart]/I");
  tree->Branch("parent", parentid, "parent[npart]/I");
  tree->Branch("stable", isstable, "isstable[npart]/I");
  tree->Branch("E", E, "E[npart]/D");
  tree->Branch("m", M, "M[npart]/D");
  tree->Branch("charge", charge, "charge[npart]/D");

  /*h_npart_outgoing = new TH1D("h_npart_outgoing", "Number of outgoing particles (proton remnants)", 100, 0., 100.);
  h_npart_outgoing_charged = new TH1D("h_npart_outgoing_charged", "Number of charged outgoing particles (proton remnants)", 100, 0., 100.);
  h_npart_outgoing_neutral = new TH1D("h_npart_outgoing_neutral", "Number of neutral outgoing particles (proton remnants)", 100, 0., 100.);*/

  int n_outgoing, n_outgoing_charged, n_outgoing_neutral;

  for (int i=0; i<100000; i++) {
    ev = *mg.GenerateOneEvent();
    //particles = ev.GetByRole(3);
    //particles = ev.GetStableParticles();
    //#ifdef DEBUG
    particles = ev.GetParticles();
    //std::cout << "--> " << particles.size() << std::endl;
    /*n_outgoing = 0;
    n_outgoing_charged = 0;
    n_outgoing_neutral = 0;*/
    for (p=particles.begin(), np=0; p!=particles.end(); p++) {
      /*if ((*p)->status==1) {
	n_outgoing++;
	if ((*p)->charge!=0.) n_outgoing_charged++;
	else n_outgoing_neutral++;
	}*/
      eta[np] = (*p)->Eta();
      phi[np] = (*p)->Phi();
      rapidity[np] = (*p)->Rapidity();
      px[np] = (*p)->px;
      py[np] = (*p)->py;
      pz[np] = (*p)->pz;
      pt[np] = (*p)->Pt();
      E[np] = (*p)->E();
      M[np] = (*p)->M();
      PID[np] = (*p)->pdgId;
      parentid[np] = (*p)->GetMother();
      isstable[np] = ((*p)->status==0 or (*p)->status==1);
      charge[np] = (*p)->charge;
      role[np] = (*p)->role;
      
      //std::cout << ":: " << (*p)->role << "\t" << (*p)->status << "\t" << (*p)->id << "\t" << (*p)->pdgId << "\t" << (*p)->name << std::endl;
      np++;
    }
    /*h_npart_outgoing->Fill(n_outgoing);
    h_npart_outgoing_charged->Fill(n_outgoing_charged);
    h_npart_outgoing_neutral->Fill(n_outgoing_neutral);*/

    tree->Fill();
    //#endif
    //ev.Dump();
    //std::cout << ev.GetLHERecord();
  }

  //tree->SaveAs("events_lpairpp_pythia.root");
  tree->SaveAs("events_lpairpp_jetset.root");

  /*gStyle->SetOptStat(0);

  c_npart = new TCanvas();
  leg = new TLegend(.65, .7, .85, .85);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);

  h_npart_outgoing_charged->Draw();
  h_npart_outgoing_charged->SetLineColor(kRed);
  leg->AddEntry(h_npart_outgoing_charged, "Charged particles");

  h_npart_outgoing_neutral->Draw("SAME");
  h_npart_outgoing_neutral->SetLineColor(kGreen);
  leg->AddEntry(h_npart_outgoing_neutral, "Neutral particles");

  h_npart_outgoing->Draw("SAME");
  leg->AddEntry(h_npart_outgoing, "All particles");

  leg->Draw();*/

  //c_npart->SaveAs("num_particles.png");

  return 0;
}

