#include "CepGen/Core/Exception.h"

using namespace cepgen;

int main()
{
  //--- try with a bit of unicode too
  const std::string test_string = "Haha, ceci est un test à géométrie variable! ☺";
  try {
    throw LoggedException( "Test", Exception::Type::warning ) << test_string;
  } catch ( const Exception& e ) {
    if ( e.message() == test_string )
      std::cout << "Test passed!" << std::endl;
      return 0;
  }
  return -1;
}
