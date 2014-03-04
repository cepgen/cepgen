#include <iostream>

#include "include/mcgen.h"

// ROOT includes
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
  TTree *tree;
  Pythia6Hadroniser had;
  //Jetset7Hadroniser had;

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
  //const int ngen = 1e4;
  const int ngen = 1e1;

  int np;
  double xsect, errxsect;
  double mx_p1, mx_p2;
  double eta[maxpart], phi[maxpart], rapidity[maxpart];
  double px[maxpart], py[maxpart], pz[maxpart], pt[maxpart], E[maxpart], M[maxpart], charge[maxpart];
  int PID[maxpart], parentid[maxpart], isstable[maxpart], role[maxpart];
  float gen_time;

  tree = new TTree("h4444", "A TTree containing information from the events produced from LPAIR++");
  tree->Branch("xsect", &xsect, "xsect/D");
  tree->Branch("errxsect", &errxsect, "errxsect/D");
  tree->Branch("MX1", &mx_p1, "MX1/D");
  tree->Branch("MX2", &mx_p2, "MX2/D");
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
  tree->Branch("gen_time", &gen_time, "gen_time/F");

  xsect = xsec;
  errxsect = err;
  for (int i=0; i<ngen; i++) {
    ev = *mg.GenerateOneEvent();
    if (i%10000==0) std::cout << "event " << i << " generated" << std::endl;
    //ev.Dump();
    //particles = ev.GetStableParticles();
    particles = ev.GetParticles();
    mx_p1 = ev.GetOneByRole(3)->M();
    mx_p2 = ev.GetOneByRole(5)->M();
    gen_time = ev.time_cpu;
    for (p=particles.begin(), np=0; p!=particles.end(); p++) {
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

    tree->Fill();
    //std::cout << ev.GetLHERecord();
  }

  tree->SaveAs("events_lpairpp_pythia.root");
  //tree->SaveAs("events_lpairpp_jetset.root");

  return 0;
}

