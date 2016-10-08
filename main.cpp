#include <iostream>

#include "core/MCGen.h"

using namespace std;

/**
 * Main caller for this Monte Carlo generator. Loads the configuration files'
 * variables if set as an argument to this program, else loads a default
 * "LHC-like" configuration, then launches the cross-section computation and
 * the events generation.
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main( int argc, char* argv[] ) {
  MCGen mg;
  
  //Logger::GetInstance()->Level = Logger::Debug;
  //Logger::GetInstance()->Level = Logger::DebugInsideLoop;
  //Logger::GetInstance()->OutputStream = ofstream( "log.txt" );
  
  if ( argc==1 ) {
    Information( "No config file provided. Setting the default parameters." );
    
    mg.parameters->process = new GamGamLL;
    mg.parameters->process_mode = Kinematics::InelasticElastic;
    mg.parameters->remnant_mode = GenericProcess::SuriYennie;

#ifdef PYTHIA6
    mg.parameters->hadroniser = new Pythia6Hadroniser;
#elifdef JETSET
    mg.parameters->hadroniser = new Jetset7Hadroniser;
#endif
    
    mg.parameters->in1p = 4000.;
    mg.parameters->in2p = 4000.;
    mg.parameters->pair = Particle::Muon;
    mg.parameters->mcut = Kinematics::BothLeptons;
    mg.parameters->minenergy = 0.; //FIXME
    mg.parameters->minpt = 5.;
    mg.parameters->mineta = -2.5;
    mg.parameters->maxeta = 2.5;
    mg.parameters->ncvg = 5e4; //FIXME
    mg.parameters->generation = true;
    //mg.parameters->maxgen = 2;
    mg.parameters->maxgen = 2e4;
  }
  else {
    Debugging( Form( "Reading config file stored in %s", argv[1] ) );
    if ( !mg.parameters->ReadConfigFile( argv[1] ) ) {
      Information( Form( "Error reading the configuration!\n\t"
                         "Please check your input file (%s)", argv[1] ) );
      return -1;
    }
  }

  // We might want to cross-check visually the validity of our run
  mg.parameters->Dump();

  // Let there be cross-section...
  double xsec, err;
  mg.ComputeXsection( &xsec, &err );

  if ( mg.parameters->generation ) {
    // The events generation starts here !
    Event ev;
    for ( unsigned int i=0; i<mg.parameters->maxgen; i++ ) {
      if ( i%1000==0 ) Information( Form( "Generating event #%d", i ) );
      ev = *mg.GenerateOneEvent();
      //ev.Dump();
    }
  }

  //mg.parameters->StoreConfigFile( "lastrun.card" );

  return 0;
}

