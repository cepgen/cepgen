#include "CepGen/Core/Exception.h"

using namespace cepgen;

int main()
{
  //--- try with a bit of unicode too
  const std::string test_string = "Haha, ceci est un test à géométrie variable! ☺";
  for ( const auto& type : { Exception::Type::debug, Exception::Type::warning, Exception::Type::error, Exception::Type::fatal } ) {
    try {
      throw LoggedException( "Test", type ) << test_string;
    } catch ( const Exception& e ) {
      if ( e.message() == test_string )
        std::cout << "Test passed for type " << (int)type << "!" << std::endl;
      e.dump();
    }
  }
  return -1;
}
