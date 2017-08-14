#include <iostream>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "TreeEvent.h"

#include "CepGen/Generator.h"
#include "CepGen/Cards/Handler.h"

using namespace std;

/**
 * Generation of events and storage in a ROOT format
 * @author Laurent Forthomme <laurent.forthomme@cern.ch>
 * @date 27 jan 2014
 */
int main( int argc, char* argv[] ) {
  const int ngen = 1e5;
  //const int ngen = 1e4;

  CepGen::Generator mg;

  double xsec, err;
  
  TFile *file;
  TTree *tree;

  CepGen::TreeEvent ev;
  
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
  ev.create( tree );
  
  ev.xsect = xsec;
  ev.errxsect = err;
  ev.litigious_events = 0;
  for ( unsigned int i=0; i<mg.parameters->maxgen; i++ ) {
    auto event = *mg.generateOneEvent();
    if ( i%10000==0 ) {
      cout << ">> event " << i << " generated" << endl;
      event.dump();
    }
    CepGen::ParticlesRef particles = event.particles();
    ev.mx_p1 = event.getOneByRole( CepGen::Particle::OutgoingBeam1 )->mass();
    ev.mx_p2 = event.getOneByRole( CepGen::Particle::OutgoingBeam2 )->mass();
    ev.hadr_trials = event.num_hadronisation_trials;

    ev.gen_time = event.time_generation;
    ev.tot_time = event.time_total;
    ev.np = 0;
    for ( const auto& p : particles ) {
      const CepGen::Particle::Momentum m = p->momentum();

      ev.kinematics[ev.np].SetXYZM( m.px(), m.py(), m.pz(), m.mass() );
      ev.rapidity[ev.np] = m.rapidity();
      ev.pt[ev.np] = m.pt();
      ev.eta[ev.np] = m.eta();
      ev.phi[ev.np] = m.phi();
      ev.E[ev.np] = p->energy();
      ev.M[ev.np] = p->mass();
      ev.PID[ev.np] = p->integerPdgId();
      ev.parentid[ev.np] = *p->mothersIds().begin();
      ev.status[ev.np] = p->status;
      ev.isstable[ev.np] = ( p->status==CepGen::Particle::Undefined or p->status==CepGen::Particle::FinalState );
      ev.charge[ev.np] = p->charge;
      ev.role[ev.np] = p->role;

      ev.np++;
    }

    tree->Fill();
  }
  cout << "Number of litigious events = " << ev.litigious_events << " -> fraction = " << ( ev.litigious_events*100./ngen ) << "%" << endl;

  file->Write();
  file->Close();
  
  return 0;
}

