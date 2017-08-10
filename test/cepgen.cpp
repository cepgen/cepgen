#include <iostream>

#include "CepGen/Generator.h"
#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Logger.h"

using namespace std;

/**
 * Main caller for this Monte Carlo generator. Loads the configuration files'
 * variables if set as an argument to this program, else loads a default
 * "LHC-like" configuration, then launches the cross-section computation and
 * the events generation.
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main( int argc, char* argv[] ) {
  CepGen::Generator mg;
  
  //CepGen::Logger::get().level = CepGen::Logger::Debug;
  //CepGen::Logger::get().level = CepGen::Logger::DebugInsideLoop;
  //CepGen::Logger::get().outputStream( ofstream( "log.txt" ) );
  
  if ( argc==1 ) {
    Information( "No config file provided. Setting the default parameters." );
    
    mg.parameters->setProcess( new CepGen::Process::GamGamLL );
    //mg.parameters->process_mode = Kinematics::InelasticElastic;
    mg.parameters->process_mode = CepGen::Kinematics::ElasticElastic;
    mg.parameters->remnant_mode = CepGen::SuriYennie;

#ifdef PYTHIA6
    mg.parameters->setHadroniser( new CepGen::Process::Pythia6Hadroniser );
#else
#ifdef JETSET
    mg.parameters->setHadroniser( new CepGen::Jetset7Hadroniser );
#endif
#endif
    
    mg.parameters->in1p = 6500.;
    mg.parameters->in2p = 6500.;
    mg.parameters->pair = CepGen::Particle::Muon;
    mg.parameters->mcut = CepGen::Kinematics::BothParticles;
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
    Information( Form( "Reading config file stored in %s", argv[1] ) );
    //CepGen::Cards::LpairReader card( argv[1] );
    CepGen::Cards::Handler<CepGen::Cards::Lpair> card( argv[1] );
    mg.setParameters( card.parameters() );
  }

  // We might want to cross-check visually the validity of our run
  mg.parameters->dump();

  // Let there be cross-section...
  double xsec, err;
  mg.computeXsection( xsec, err );

  if ( mg.parameters->generation ) {
    // The events generation starts here !
    CepGen::Event ev;
    for ( unsigned int i=0; i<mg.parameters->maxgen; i++ ) {
      ev = *mg.generateOneEvent();
      if ( i%1000==0 ) {
        Information( Form( "Generating event #%d", i ) );
        ev.dump();
      }
    }
  }

  //mg.parameters->storeConfigFile( "lastrun.card" );

  return 0;
}

