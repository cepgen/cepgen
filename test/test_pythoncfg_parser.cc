#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Parameters.h"

#include <string>
#include <iostream>

#include "ArgumentsParser.h"

using namespace std;
using namespace cepgen;

int
main( int argc, char* argv[] )
{
  string card;

  ArgumentsParser( argc, argv )
    .addArgument( "card", "input card", &card, 'i' )
    .parse();

  try {
    CG_INFO( "main" )
      << &card::PythonHandler( card.c_str() ).parameters();
  } catch ( const Exception& e ) {
    e.dump();
  }
  return 0;
}
