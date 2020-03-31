#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Cards/Handler.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  namespace card
  {
    Parameters&
    Handler::parse( const std::string& filename )
    {
      try {
        Parameters params;
        auto parser = CardsHandlerFactory::get().build( extension( filename ) );
        return parser->parse( filename, params );
      } catch ( const std::invalid_argument& err ) {
        throw CG_FATAL( "Cards:handler" )
          << "Failed to parse the steering card at \"" << filename << "\"!\n"
          << err.what();
      }
    }
  }
}
