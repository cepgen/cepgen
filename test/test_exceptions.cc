#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Logger.h"

using namespace cepgen;

int main()
{
  utils::Logger::get().level = utils::Logger::Level::nothing;
  utils::Logger::get().output = nullptr;
  //--- try with a bit of unicode too
  const std::string test_string = "Haha, ceci est un test à géométrie variable! ☺";
  for ( int type = (int)Exception::Type::undefined; type <= (int)Exception::Type::fatal; ++type ) {
    try {
      throw LoggedException( "Test", (Exception::Type)type ) << test_string;
      return -1;
    } catch ( const Exception& e ) {
      if ( e.message() == test_string )
        std::cout << "Test passed for type " << type << "!" << std::endl;
    }
  }
  return 0;
}
