#include <iostream>


// ROOT includes
#include "TFile.h"
#include "TTree.h"

#include "include/MCGen.h"

using namespace std;

/**
 * Generation of events and storage in a ROOT format
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * @date 27 jan 2014
 */
int main(int argc, char* argv[]) {
  const Int_t maxpart = 500;

  //const int ngen = 1e5;
  const int ngen = 1e4;

  MCGen mg;
  Event ev;
  //PPtoLL proc;
  //Jetset7Hadroniser had;

  double xsec, err;
  ParticlesRef particles, remn;
  ParticlesRef::iterator p;
  
  TFile *file;
  TTree *tree;
  
  int np;
  double xsect, errxsect;
  double mx_p1, mx_p2;
  double eta[maxpart], phi[maxpart], rapidity[maxpart];
  double px[maxpart], py[maxpart], pz[maxpart], pt[maxpart], E[maxpart], M[maxpart], charge[maxpart];
  int PID[maxpart], parentid[maxpart], isstable[maxpart], role[maxpart], status[maxpart];
  float gen_time, tot_time;
  int nremn_ch[2], nremn_nt[2];
  int hadr_trials, litigious_events;
 
  TString filename = "events.root";
  if (argc>3) filename = TString(argv[2]);
  file = new TFile(filename, "recreate");
  if (!file) {
    cout << "ERROR while trying to create the output file!" << endl;
  }
  if (argc==1) {
    mg.parameters->process = new PPtoLL;
    mg.parameters->in1p = 3500.;
    mg.parameters->in2p = 3500.;
    mg.parameters->pair = Particle::Muon;
    mg.parameters->mcut = 2;
    mg.parameters->minenergy = 0.; //FIXME
    mg.parameters->minpt = 5.;
    mg.parameters->maxgen = ngen;
    mg.parameters->hadroniser = new Pythia6Hadroniser;
    mg.parameters->remnant_mode = GenericProcess::SuriYennie;
    mg.parameters->process_mode = GenericProcess::ElasticElastic;
    //mg.parameters->ncvg = 5e3; //FIXME
    //mg.parameters->SetEtaRange(-2.5, 2.5);
    //mg.parameters->SetEtaRange(-999., 999.);
  }
  else {
    Debug(Form("Reading config file stored in %s", argv[1]));
    if (!mg.parameters->ReadConfigFile(argv[1])) {
      Info(Form("Error reading the configuration!\n\t"
                "Please check your input file (%s)", argv[1]));
      return -1;
    }
  }
    
  mg.parameters->generation = true;
  mg.parameters->Dump();

  mg.ComputeXsection(&xsec, &err);

  tree = new TTree("h4444", "A TTree containing information from the events produced from LPAIR++");
  tree->Branch("xsect", &xsect, "xsect/D");
  tree->Branch("errxsect", &errxsect, "errxsect/D");
  tree->Branch("MX1", &mx_p1, "MX1/D");
  tree->Branch("MX2", &mx_p2, "MX2/D");
  tree->Branch("ip", &np, "npart/I");
  tree->Branch("nremn_charged", nremn_ch, "nremn_charged[2]/I");
  tree->Branch("nremn_neutral", nremn_nt, "nremn_neutral[2]/I");
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
  tree->Branch("status", status, "status[npart]/I");
  tree->Branch("stable", isstable, "isstable[npart]/I");
  tree->Branch("E", E, "E[npart]/D");
  tree->Branch("m", M, "M[npart]/D");
  tree->Branch("charge", charge, "charge[npart]/D");
  tree->Branch("generation_time", &gen_time, "gen_time/F");
  tree->Branch("total_time", &tot_time, "gen_time/F");
  tree->Branch("hadronisation_trials", &hadr_trials, "hadronisation_trials/I");
  //tree->Branch("event", &ev);

  xsect = xsec;
  errxsect = err;
  litigious_events = 0;
  for (int i=0; i<mg.parameters->maxgen; i++) {
    ev = *mg.GenerateOneEvent();
    if (i%10000==0) {
      cout << ">> event " << i << " generated" << endl;
      ev.Dump();
    }
    //particles = ev.GetStableParticles();
    particles = ev.GetParticles();
    mx_p1 = ev.GetOneByRole(Particle::OutgoingBeam1)->M();
    mx_p2 = ev.GetOneByRole(Particle::OutgoingBeam2)->M();
    hadr_trials = ev.num_hadronisation_trials;

    /*// study the remnants
    for (unsigned int j=0; j<2; j++) { nremn_ch[j] = nremn_nt[j] = 0; }
    
    // Proton 1
    remn = ev.GetByRole(Particle::OutgoingBeam1);
    for (p=remn.begin(); p!=remn.end(); p++) {
      const Particle::Status s = (*p)->status;
      if ( s!=Particle::Undefined
       and s!=Particle::FinalState) continue; // we only count stable particles
      if ((int)(*p)->charge!=(*p)->charge) continue; // to get rid of partons
      if (s==Particle::Undecayed or s==Particle::Incoming) continue; // to get rid of partons
      
      // split between charged and neutral remnants
      if ((int)(*p)->charge%2!=0) nremn_ch[0]++;
      else nremn_nt[0]++;
    }

    // Proton 2
    remn = ev.GetByRole(Particle::OutgoingBeam2);
    for (p=remn.begin(); p!=remn.end(); p++) {
      const Particle::Status s = (*p)->status;
      if ( s!=Particle::Undefined
       and s!=Particle::FinalState) continue; // we only count stable particles
      if ((int)(*p)->charge!=(*p)->charge) continue; // to get rid of partons
      if (s==Particle::Undecayed or s==Particle::Incoming) continue; // to get rid of partons
      
      // split between charged and neutral remnants
      if ((int)(*p)->charge%2!=0) nremn_ch[1]++;
      else nremn_nt[1]++;
    }

    if (nremn_ch[0]%2==0 or nremn_ch[1]%2==0) {
      //ev.Dump();
      cout << "--> Event " << i << " contains" << endl
           << "\t-> Remnants 1: " << nremn_ch[0] << " charged and " << nremn_nt[0] << " neutral remnants" << endl
           << "\t-> Remnants 2: " << nremn_ch[1] << " charged and " << nremn_nt[1] << " neutral remnants" << endl;
      litigious_events++;
    }*/
    gen_time = ev.time_generation;
    tot_time = ev.time_total;
    for (p=particles.begin(), np=0; p!=particles.end(); p++) {
      const Particle::Momentum m = (*p)->GetMomentum();
      eta[np] = m.Eta();
      phi[np] = m.Phi();
      rapidity[np] = m.Rapidity();
      px[np] = m.Px();
      py[np] = m.Py();
      pz[np] = m.Pz();
      pt[np] = m.Pt();
      E[np] = (*p)->E();
      M[np] = (*p)->M();
      PID[np] = (*p)->GetIntPDGId();
      parentid[np] = *(*p)->GetMothersIds().begin();
      status[np] = (*p)->status;
      isstable[np] = ((*p)->status==Particle::Undefined or (*p)->status==Particle::FinalState);
      charge[np] = (*p)->charge;
      role[np] = (*p)->role;
      
      //cout << ":: " << (*p)->role << "\t" << (*p)->status << "\t" << (*p)->id << "\t" << (*p)->pdgId << "\t" << (*p)->name << endl;
      np++;
    }

    tree->Fill();
    //cout << ev.GetLHERecord();
  }
  cout << "Number of litigious events = " << litigious_events << " -> fraction = " << (double)litigious_events/ngen*100 << "%" << endl;

  //tree->SaveAs("events.root");
  //tree->SaveAs("events_lpairpp_pythia.root");
  //tree->SaveAs("events_lpairpp_jetset.root");
  //tree->SaveAs("events_lpairpp_elastic_pt5.root");
  //tree->SaveAs("events_lpairpp_singlediss_pythia_pt5.root");
  //tree->SaveAs("events_lpairpp_doublediss_pythia_pt5.root");
  file->Write();
  file->Close();
  
  //delete file;
  //delete tree;

  return 0;
}

