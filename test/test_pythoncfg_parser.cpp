#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Core/Exception.h"

#include <string>
#include <iostream>

using namespace std;

int
main( int argc, char* argv[] )
{
  if ( argc < 2 )
    throw CG_FATAL( "main" ) << "One argument required!";

  try {
    CepGen::cards::PythonHandler py( argv[1] );
    CG_INFO( "main" ) << py.parameters();
  } catch ( CepGen::Exception& e ) {
    e.dump();
  }
  return 0;
}
