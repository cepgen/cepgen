#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Parameters.h"

#include <string>
#include <iostream>

using namespace std;
using namespace cepgen;

int
main( int argc, char* argv[] )
{
  if ( argc < 2 )
    throw CG_FATAL( "main" ) << "One argument required!";

  try {
    CG_INFO( "main" )
      << &card::PythonHandler( argv[1] ).parameters();
  } catch ( const Exception& e ) {
    e.dump();
  }
  return 0;
}
