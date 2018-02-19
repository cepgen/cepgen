#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Core/Exception.h"

#include <string>
#include <iostream>

using namespace std;

int
main( int argc, char* argv[] )
{
  if ( argc < 2 )
    FatalError( "One argument required!" );

  CepGen::Cards::PythonHandler py( argv[1] );
  py.parameters().dump();
  return 0;
}
