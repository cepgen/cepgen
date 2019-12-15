#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Cards/Handler.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  namespace card
  {
    std::unique_ptr<Handler>
    Handler::parse( const std::string& filename )
    {
      try {
        return CardsHandlerFactory::get().build( extension( filename ),
          ParametersList().set<std::string>( FILENAME_KEY, filename ) );
      } catch ( const std::invalid_argument& err ) {
        throw CG_FATAL( "Cards:handler" )
          << "Failed to parse the steering card at \"" << filename << "\"!\n"
          << err.what();
      }
    }
  }
}
