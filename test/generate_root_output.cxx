#include <iostream>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "CepGen/Generator.h"
#include "CepGen/Cards/Handler.h"

using namespace std;

/**
 * Generation of events and storage in a ROOT format
 * @author Laurent Forthomme <laurent.forthomme@cern.ch>
 * @date 27 jan 2014
 */
int main( int argc, char* argv[] ) {
  const Int_t maxpart = 10;

  const int ngen = 1e5;
  //const int ngen = 1e4;

  CepGen::Generator mg;
  CepGen::Event ev;

  double xsec, err;
  
  TFile *file;
  TTree *tree;
  
  int np;
  double xsect, errxsect;
  double mx_p1, mx_p2;
  double pt[maxpart], eta[maxpart], phi[maxpart], rapidity[maxpart];
  double E[maxpart], M[maxpart], charge[maxpart];
  int PID[maxpart], parentid[maxpart], isstable[maxpart], role[maxpart], status[maxpart];
  TLorentzVector kinematics[maxpart];
  float gen_time, tot_time;
  int nremn_ch[2], nremn_nt[2];
  int hadr_trials, litigious_events;
 
  if ( argc<2 ) {
    InError( Form( "Usage: %s <input card> [output .root filename]\n\t"
                   "   or: %s <process type 1..4> [output .root filename]", argv[0], argv[0] ) );
    return -1;
  }

  if ( atoi( argv[1] )<=4 and atoi( argv[1] )>0 ) {
    // do not provide an input card
    mg.parameters->setProcess( new CepGen::Process::GamGamLL );
    mg.parameters->in1p = 6500.;
    mg.parameters->in2p = 6500.;
    mg.parameters->pair = CepGen::Particle::Muon;
    mg.parameters->mcut = CepGen::Kinematics::BothParticles;
    mg.parameters->minenergy = 0.; //FIXME
    mg.parameters->minpt = 15.;
    mg.parameters->maxgen = ngen;
    mg.parameters->remnant_mode = CepGen::SuriYennie;
    mg.parameters->process_mode = ( argc>1 )
      ? static_cast<CepGen::Kinematics::ProcessMode>( atoi( argv[1] ) )
      : CepGen::Kinematics::ElasticElastic;
    //mg.parameters->ncvg = 5e3; //FIXME
    mg.parameters->mineta = -2.5;
    mg.parameters->maxeta = 2.5;
    mg.parameters->maxmx = 1.e3;
  }
  else {
    mg.setParameters( CepGen::Cards::LpairReader( argv[1] ).parameters() );
  }
    
  mg.parameters->generation = true;
  mg.parameters->dump();

  //----- open the output root file

  const TString filename = ( argc>2 ) ? argv[2] : "events.root";
  file = new TFile( filename, "recreate" );
  if ( !file ) {
    cout << "ERROR while trying to create the output file!" << endl;
  }

  //----- start by computing the cross section for the list of parameters applied
  mg.computeXsection( xsec, err );

  tree = new TTree( "h4444", "A TTree containing information from the events produced from CepGen" );
  tree->Branch( "xsect", &xsect, "xsect/D" );
  tree->Branch( "errxsect", &errxsect, "errxsect/D" );
  tree->Branch( "MX1", &mx_p1, "MX1/D" );
  tree->Branch( "MX2", &mx_p2, "MX2/D" );
  tree->Branch( "npart", &np, "npart/I" );
  tree->Branch( "nremn_charged", nremn_ch, "nremn_charged[2]/I" );
  tree->Branch( "nremn_neutral", nremn_nt, "nremn_neutral[2]/I" );
  //tree->Branch( "kinematics", kinematics, "TLorentzVector[npart]" );
  tree->Branch( "rapidity", rapidity, "rapidity[npart]/D" );
  tree->Branch( "pt", pt, "pt[npart]/D" );
  tree->Branch( "eta", eta, "eta[npart]/D" );
  tree->Branch( "phi", phi, "phi[npart]/D" );
  tree->Branch( "icode", PID, "icode[npart]/I" );
  tree->Branch( "role", role, "role[npart]/I" );
  tree->Branch( "parent", parentid, "parent[npart]/I" );
  tree->Branch( "status", status, "status[npart]/I" );
  tree->Branch( "stable", isstable, "stable[npart]/I" );
  tree->Branch( "E", E, "E[npart]/D" );
  tree->Branch( "m", M, "m[npart]/D" );
  tree->Branch( "charge", charge, "charge[npart]/D" );
  tree->Branch( "generation_time", &gen_time, "generation_time/F" );
  tree->Branch( "total_time", &tot_time, "total_time/F" );
  tree->Branch( "hadronisation_trials", &hadr_trials, "hadronisation_trials/I" );
  
  xsect = xsec;
  errxsect = err;
  litigious_events = 0;
  for ( unsigned int i=0; i<mg.parameters->maxgen; i++ ) {
    ev = *mg.generateOneEvent();
    if ( i%10000==0 ) {
      cout << ">> event " << i << " generated" << endl;
      ev.dump();
    }
    CepGen::ParticlesRef particles = ev.particles();
    mx_p1 = ev.getOneByRole( CepGen::Particle::OutgoingBeam1 )->mass();
    mx_p2 = ev.getOneByRole( CepGen::Particle::OutgoingBeam2 )->mass();
    hadr_trials = ev.num_hadronisation_trials;

    gen_time = ev.time_generation;
    tot_time = ev.time_total;
    np = 0;
    for ( const auto& p : particles ) {
      const CepGen::Particle::Momentum m = p->momentum();

      kinematics[np].SetXYZM( m.px(), m.py(), m.pz(), m.mass() );
      rapidity[np] = m.rapidity();
      pt[np] = m.pt();
      eta[np] = m.eta();
      phi[np] = m.phi();
      E[np] = p->energy();
      M[np] = p->mass();
      PID[np] = p->integerPdgId();
      parentid[np] = *p->mothersIds().begin();
      status[np] = p->status;
      isstable[np] = ( p->status==CepGen::Particle::Undefined or p->status==CepGen::Particle::FinalState );
      charge[np] = p->charge;
      role[np] = p->role;
      
      np++;
    }

    tree->Fill();
  }
  cout << "Number of litigious events = " << litigious_events << " -> fraction = " << static_cast<double>( litigious_events )/ngen*100 << "%" << endl;

  file->Write();
  file->Close();
  
  return 0;
}

