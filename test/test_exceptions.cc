#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Logger.h"
#include <fstream>

using namespace cepgen;

int main()
{
  utils::Logger::get().level = utils::Logger::Level::nothing;
  //utils::Logger::get().output = new std::ofstream( "test.log" );
  utils::Logger::get().output = nullptr;

  //--- try with a bit of unicode too
  const std::string test_string = "Haha, ceci est un test à géométrie variable! ☺";
  for ( int type = (int)Exception::Type::undefined; type < (int)Exception::Type::fatal; ++type ) {
    try {
      throw LoggedException( "Test", (Exception::Type)type ) << test_string;
      std::cout << "Test failed for type " << type << "!" << std::endl;
      return -1;
    } catch ( const Exception& e ) {
      if ( e.message() == test_string )
        std::cout << "Test passed for type " << type << "!" << std::endl;
      else
        std::cout << "Test passed for type " << type << " (unicode)!" << std::endl;
    }
  }
  return 0;
}
