#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Cards/LpairHandler.h"

#include "CepGen/Generator.h"

#include "CepGen/Core/Logger.h"
#include "CepGen/Core/Exception.h"

// necessary includes to build the default run
#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <iostream>

using namespace std;

void printEvent( const CepGen::Event& ev, unsigned long ev_id )
{
//cout << ev_id << endl;
  if ( ev_id % 5000 != 0 )
    return;

  Information( Form( "Generating event #%d", ev_id ) );
  ev.dump();
}

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
    mg.parameters->kinematics.inp = { 6500., 6500. };
    mg.parameters->kinematics.central_system = { CepGen::Muon, CepGen::Muon };
    mg.parameters->kinematics.cuts.central[CepGen::Cuts::pt_single].min() = 15.;
    mg.parameters->kinematics.cuts.central[CepGen::Cuts::eta_single] = { -2.5, 2.5 };
    mg.parameters->integrator.ncvg = 5e4;
    mg.parameters->generation.num_threads = 4;
    mg.parameters->generation.enabled = true;
    mg.parameters->generation.maxgen = 1e5;
  }
  else {
    Information( Form( "Reading config file stored in %s", argv[1] ) );
    //CepGen::Cards::LpairReader card( argv[1] );
    const std::string extension = CepGen::Cards::Handler::getExtension( argv[1] );
    if ( extension == "card" )
      mg.setParameters( CepGen::Cards::LpairHandler( argv[1] ).parameters() );
#ifdef PYTHON
    else if ( extension == "py" )
      mg.setParameters( CepGen::Cards::PythonHandler( argv[1] ).parameters() );
#endif
  }

  // We might want to cross-check visually the validity of our run
  mg.parameters->dump();

  // Let there be cross-section...
  double xsec = 0., err = 0.;
  mg.computeXsection( xsec, err );

  if ( mg.parameters->generation.enabled )
    // The events generation starts here !
    mg.generate( printEvent ); // use a callback function!

  return 0;
}

