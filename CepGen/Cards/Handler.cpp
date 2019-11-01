#include "CepGen/Cards/CardsHandler.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  namespace card
  {
    Parameters&
    Handler::parse( const char* filename )
    {
      try {
        static auto hnd = CardsHandler::get().build(
          extension( filename ),
          ParametersList().set<std::string>( FILENAME_KEY, filename ) );
        return hnd->parameters();
      } catch ( const std::invalid_argument& err ) {
        throw CG_FATAL( "Cards:handler" )
          << "Failed to parse the steering card at \"" << filename << "\"!\n"
          << err.what();
      }
    }
  }
}
