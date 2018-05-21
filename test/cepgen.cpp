#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Cards/LpairHandler.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Generator.h"

#include "CepGen/Core/Logger.h"
#include "CepGen/Core/Exception.h"

// necessary includes to build the default run
#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <iostream>
#include <fstream>

using namespace std;

void printEvent( const CepGen::Event& ev, unsigned long ev_id )
{
  if ( ev_id % 5000 != 0 )
    return;

  CG_INFO( "printEvent" ) << "Generating event #" << ev_id << ".";
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
  //CepGen::Logger::get( new ofstream( "log.txt" ) );
  //CepGen::Logger::get().level = CepGen::Logger::Level::Debug;

  CepGen::Generator mg;

  if ( argc > 1 ) {
    CG_INFO( "main" ) << "Reading config file stored in " << argv[1] << ".";
    const std::string extension = CepGen::Cards::Handler::getExtension( argv[1] );
    if ( extension == "card" )
      mg.setParameters( CepGen::Cards::LpairHandler( argv[1] ).parameters() );
#ifdef PYTHON
    else if ( extension == "py" )
      mg.setParameters( CepGen::Cards::PythonHandler( argv[1] ).parameters() );
#endif
    else
      throw CG_FATAL( "main" ) << "Unrecognized steering card extension: ." << extension << "!";
  }
  else {
    CG_INFO( "main" ) << "No config file provided. Setting the default parameters.";

    mg.parameters->setProcess( new CepGen::Process::GamGamLL );
    mg.parameters->kinematics.mode = CepGen::Kinematics::Mode::ElasticElastic;
    mg.parameters->kinematics.structure_functions = CepGen::StructureFunctions::SuriYennie;
    mg.parameters->kinematics.inp = { 6500., 6500. };
    mg.parameters->kinematics.central_system = { CepGen::PDG::Muon, CepGen::PDG::Muon };
    mg.parameters->kinematics.cuts.central.pt_single.min() = 15.;
    mg.parameters->kinematics.cuts.central.eta_single = { -2.5, 2.5 };
    mg.parameters->integrator.ncvg = 5e4;
    mg.parameters->generation.num_threads = 4;
    mg.parameters->generation.enabled = true;
    mg.parameters->generation.maxgen = 1e3;
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

