#include <iostream>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "TreeEvent.h"

#include "CepGen/Generator.h"
#include "CepGen/Cards/LpairHandler.h"
#include "CepGen/Cards/ConfigHandler.h"

using namespace std;

/**
 * Generation of events and storage in a ROOT format
 * @author Laurent Forthomme <laurent.forthomme@cern.ch>
 * @date 27 jan 2014
 */
int main( int argc, char* argv[] ) {
  CepGen::Generator mg;

  if ( argc < 2 ) {
    InError( Form( "Usage: %s <input card> [output .root filename]", argv[0] ) );
    return -1;
  }
  const std::string incard( argv[1] ), extension = incard.substr( incard.find_last_of( "." )+1 );
  if ( extension == "card" ) mg.setParameters( CepGen::Cards::LpairHandler( argv[1] ).parameters() );
  else if ( extension == "cfg" ) mg.setParameters( CepGen::Cards::ConfigHandler( argv[1] ).parameters() );

  mg.parameters->generation.enabled = true;
  mg.parameters->dump();

  //----- open the output root file

  const TString filename = ( argc > 2 ) ? argv[2] : "events.root";
  auto file = TFile::Open( filename, "recreate" );
  if ( !file ) {
    cout << "ERROR while trying to create the output file!" << endl;
  }

  //----- start by computing the cross section for the list of parameters applied
  double xsec, err;
  mg.computeXsection( xsec, err );

  //----- then generate the events and the container tree structure

  auto run_tree = new TTree( "run", "A TTree containing information from the run produced with CepGen" );
  auto ev_tree = new TTree( "events", "A TTree containing information from the events produced from CepGen" );

  CepGen::TreeRun run;
  run.create( run_tree );
  run.xsect = xsec;
  run.errxsect = err;
  run.litigious_events = 0;
  run_tree->Fill();

  CepGen::TreeEvent ev;
  ev.create( ev_tree );

  for ( unsigned int i = 0; i < mg.parameters->generation.maxgen; ++i ) {
    const auto event = mg.generateOneEvent();
    if ( !event ) FatalError( "Failed to generate the event!" );

    ev.clear();
    if ( i % 10000 == 0 ) {
      cout << ">> event " << i << " generated" << endl;
      //event->dump();
    }

    ev.gen_time = event->time_generation;
    ev.tot_time = event->time_total;
    ev.np = 0;
    for ( const auto& p : event->particles() ) {
      const CepGen::Particle::Momentum m = p.momentum();

      //ev.kinematics[ev.np].SetXYZM( m.px(), m.py(), m.pz(), m.mass() );
      ev.rapidity[ev.np] = m.rapidity();
      ev.pt[ev.np] = m.pt();
      ev.eta[ev.np] = m.eta();
      ev.phi[ev.np] = m.phi();
      ev.E[ev.np] = p.energy();
      ev.m[ev.np] = p.mass();
      ev.pdg_id[ev.np] = p.integerPdgId();
      ev.parent1[ev.np] = ( p.mothers().size() > 0 ) ? *p.mothers().begin() : -1;
      ev.parent2[ev.np] = ( p.mothers().size() > 1 ) ? *p.mothers().rbegin() : -1;
      ev.status[ev.np] = p.status();
      ev.stable[ev.np] = ( (short)p.status() > 0 );
      ev.charge[ev.np] = p.charge();
      ev.role[ev.np] = p.role();

      ev.np++;
    }

    ev_tree->Fill();
  }
  //cout << "Number of litigious events = " << run.litigious_events << " -> fraction = " << ( run.litigious_events*100./ngen ) << "%" << endl;

  file->Write();
  delete file;

  return 0;
}

