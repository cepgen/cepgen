//--- steering cards
#include "CepGen/Cards/Handler.h"

#include "CepGen/Generator.h"
#include "CepGen/Core/Exception.h"

#include "AbortHandler.h"

using namespace std;

/** Example executable for CepGen
 * - loads the steering card variables into the environment,
 * - launches the cross-section computation and the events generation (if requested).
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \addtogroup Executables
 */
int main( int argc, char* argv[] )
{
  if ( argc < 2 )
    throw CG_FATAL( "main" )
      << "No config file provided!\n\t"
      << "Usage: " << argv[0] << " config-file";

  //--- first start by defining the generator object
  cepgen::Generator gen;
  gen.setParameters( cepgen::card::Handler::parse( argv[1] ) );

  //--- list all parameters
  CG_LOG( "main" ) << gen.parametersPtr();

  cepgen::utils::AbortHandler ctrl_c;

  try {
    //--- let there be a cross-section...
    double xsec = 0., err = 0.;
    gen.computeXsection( xsec, err );

    if ( gen.parameters().generation().enabled )
      //--- events generation starts here
      // (one may use a callback function)
      gen.generate();
  } catch ( const cepgen::utils::RunAbortedException& e ) {
    CG_DEBUG( "main" ) << "Run aborted!";
  } catch ( const cepgen::Exception& e ) {
    e.dump();
  } catch ( const std::exception& e ) {
    CG_FATAL( "main" ) << "Other exception caught!\n\t"
      << e.what();
  }

  return 0;
}
