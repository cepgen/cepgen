#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Cards/LpairHandler.h"
#include "CepGen/Generator.h"
#include "CepGen/Event/Event.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "TreeInfo.h"
#include "abort.h"

#include <iostream>

using namespace std;

std::unique_ptr<CepGen::TreeRun> run;
std::unique_ptr<CepGen::TreeEvent> ev;

void fill_event_tree( const CepGen::Event& event, unsigned long ev_id )
{
  //if ( ev_id % 10 == 0 )
  //  cout << ">> event " << ev_id << " generated" << endl;

  if ( !ev || !run )
    return;

  ev->gen_time = event.time_generation;
  ev->tot_time = event.time_total;
  ev->np = 0;
  //cout << event.particles().size() << endl;
  ev->momentum.reserve( event.particles().size() );
  for ( const auto& p : event.particles() ) {
    const CepGen::Particle::Momentum m = p.momentum();
    ev->momentum[ev->np].SetPxPyPzE( m.px(), m.py(), m.pz(), p.energy() );
    ev->rapidity[ev->np] = m.rapidity();
    ev->pt[ev->np] = m.pt();
    ev->eta[ev->np] = m.eta();
    ev->phi[ev->np] = m.phi();
    ev->E[ev->np] = p.energy();
    ev->m[ev->np] = p.mass();
    ev->pdg_id[ev->np] = p.integerPdgId();
    ev->parent1[ev->np] = ( p.mothers().size() > 0 ) ? *p.mothers().begin() : -1;
    ev->parent2[ev->np] = ( p.mothers().size() > 1 ) ? *p.mothers().rbegin() : -1;
    ev->status[ev->np] = (int)p.status();
    ev->stable[ev->np] = ( (short)p.status() > 0 );
    ev->charge[ev->np] = p.charge();
    ev->role[ev->np] = p.role();

    ev->np++;
  }
  run->num_events += 1;
  ev->fill();
}

/**
 * Generation of events and storage in a ROOT format
 * @author Laurent Forthomme <laurent.forthomme@cern.ch>
 * @date 27 jan 2014
 */
int main( int argc, char* argv[] ) {
  CepGen::Generator mg;

  if ( argc < 2 )
    throw CG_FATAL( "main" ) << "Usage: " << argv[0] << " input-card [filename=events.root]";

  const std::string extension = CepGen::Cards::Handler::getExtension( argv[1] );
  if ( extension == "card" )
    mg.setParameters( CepGen::Cards::LpairHandler( argv[1] ).parameters() );
  else if ( extension == "py" )
    mg.setParameters( CepGen::Cards::PythonHandler( argv[1] ).parameters() );

  mg.parameters->generation.enabled = true;
  mg.parameters->dump();

  //----- open the output root file

  const char* filename = ( argc > 2 ) ? argv[2] : "events.root";
  std::unique_ptr<TFile> file( TFile::Open( filename, "recreate" ) );
  if ( !file )
    throw CG_FATAL( "main" ) << "Failed to create the output file!";

  AbortHandler ctrl_c;
  //----- start by computing the cross section for the list of parameters applied
  double xsec, err;
  mg.computeXsection( xsec, err );

  //----- then generate the events and the container tree structure

  std::unique_ptr<TTree> ev_tree( new TTree( "events", "A TTree containing information from the events produced from CepGen" ) );

  run.reset( new CepGen::TreeRun );
  run->create();
  run->xsect = xsec;
  run->errxsect = err;
  run->litigious_events = 0;
  run->sqrt_s = mg.parameters->kinematics.sqrtS();
  run->fill();

  ev.reset( new CepGen::TreeEvent );
  ev->create( ev_tree.get() );

  // launch the events generation
  try {
    mg.generate( fill_event_tree );
  } catch ( const CepGen::Exception& ) {}

  file->Write();
  CG_INFO( "main" ) << "Events written on \"" << filename << "\".";

  return 0;
}
