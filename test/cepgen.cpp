#include <iostream>

#include "CepGen/Generator.h"

#include "CepGen/Cards/LpairHandler.h"
#include "CepGen/Cards/ConfigHandler.h"
#include "CepGen/Core/Logger.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Hadronisers/Pythia8Hadroniser.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"

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
  
  if ( argc == 1 ) {
    Information( "No config file provided. Setting the default parameters." );
    
    mg.parameters->setProcess( new CepGen::Process::GamGamLL );
    //mg.parameters->process_mode = Kinematics::InelasticElastic;
    mg.parameters->kinematics.mode = CepGen::Kinematics::ElasticElastic;
    mg.parameters->kinematics.structure_functions = CepGen::StructureFunctions::SuriYennie;
#ifdef PYTHIA8
    mg.parameters->setHadroniser( new CepGen::Hadroniser::Pythia8Hadroniser );
#endif

    mg.parameters->kinematics.inp = { 6500., 6500. };
    mg.parameters->kinematics.central_system = { CepGen::Muon, CepGen::Muon };
    mg.parameters->kinematics.cuts.central[CepGen::Cuts::pt_single].min() = 15.;
    mg.parameters->kinematics.cuts.central[CepGen::Cuts::eta_single] = { -2.5, 2.5 };
    mg.parameters->integrator.ncvg = 5e4; //FIXME
    mg.parameters->generation.enabled = true;
    //mg.parameters->maxgen = 2;
    mg.parameters->generation.maxgen = 2e4;
  }
  else {
    Information( Form( "Reading config file stored in %s", argv[1] ) );
    //CepGen::Cards::LpairReader card( argv[1] );
    const std::string file( argv[1] ), extension = file.substr( file.find_last_of( "." )+1 );
    if ( extension == "card" ) mg.setParameters( CepGen::Cards::LpairHandler( argv[1] ).parameters() );
    else if ( extension == "cfg" ) mg.setParameters( CepGen::Cards::ConfigHandler( argv[1] ).parameters() );
  }

  // We might want to cross-check visually the validity of our run
  mg.parameters->dump();

  // Let there be cross-section...
  double xsec, err;
  mg.computeXsection( xsec, err );

  if ( mg.parameters->generation.enabled ) {
    // The events generation starts here !
    CepGen::Event ev;
    for ( unsigned int i=0; i<mg.parameters->generation.maxgen; i++ ) {
      ev = *mg.generateOneEvent();
      if ( i%1000==0 ) {
        Information( Form( "Generating event #%d", i ) );
        ev.dump();
      }
    }
  }

  // store the current configuration
  CepGen::Cards::ConfigHandler::store( mg.parameters.get(), "last_run.cfg" );

  return 0;
}

