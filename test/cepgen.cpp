//--- steering cards
#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Cards/LpairHandler.h"

#include "CepGen/Generator.h"
#include "CepGen/Core/Exception.h"

//--- necessary include to build the default run
#include "CepGen/Physics/PDG.h"
#include "CepGen/Processes/ProcessesHandler.h"

#include "AbortHandler.h"

#include <iostream>

using namespace std;

/**
 * Main caller for this MC generator.
 *  * loads the configuration files' variables if passed as an argument,
 *    or a default LPAIR-like configuration,
 *  * launches the cross-section computation and the events generation.
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main( int argc, char* argv[] )
{
  //--- first start by defining the generator object
  cepgen::Generator gen;

  if ( argc < 2 ) {
    CG_INFO( "main" ) << "No config file provided. Setting the default parameters.";

    //--- default run: LPAIR elastic ɣɣ → µ⁺µ¯ at 13 TeV
    cepgen::ParametersList pgen;
    pgen.set<int>( "pair", (int)cepgen::PDG::muon );
    gen.parameters->setProcess( cepgen::proc::ProcessesHandler::get().build( "lpair", pgen ) );
    gen.parameters->kinematics.mode = cepgen::KinematicsMode::ElasticElastic;
    gen.parameters->kinematics.cuts.central.pt_single.min() = 15.;
    gen.parameters->kinematics.cuts.central.eta_single = { -2.5, 2.5 };
    gen.parameters->generation.enabled = true;
    gen.parameters->generation.maxgen = 1e3;
  }
  else {
    CG_INFO( "main" ) << "Reading config file stored in " << argv[1] << ".";
    const std::string extension = cepgen::card::Handler::getExtension( argv[1] );
    if ( extension == "card" )
      gen.setParameters( cepgen::card::LpairHandler( argv[1] ).parameters() );
#ifdef PYTHON
    else if ( extension == "py" )
      gen.setParameters( cepgen::card::PythonHandler( argv[1] ).parameters() );
#endif
    else
      throw CG_FATAL( "main" ) << "Unrecognized steering card extension: ." << extension << "!";
  }

  //--- list all parameters
  CG_LOG( "main" ) << gen.parameters.get();

  cepgen::utils::AbortHandler ctrl_c;

  try {
    //--- let there be a cross-section...
    double xsec = 0., err = 0.;
    gen.computeXsection( xsec, err );

    if ( gen.parameters->generation.enabled )
      //--- events generation starts here
      // (one may use a callback function)
      gen.generate();
  } catch ( const cepgen::utils::RunAbortedException& e ) {
    CG_INFO( "main" ) << "Run aborted!";
  } catch ( const cepgen::Exception& e ) {
    e.dump();
  }

  return 0;
}
