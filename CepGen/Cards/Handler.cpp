#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Cards/Handler.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  namespace card
  {
    Parameters*
    Handler::parse( const std::string& filename )
    {
      try {
        auto parser = CardsHandlerFactory::get().build( extension( filename ) );
        return parser->parse( filename, new Parameters );
      } catch ( const std::invalid_argument& err ) {
        throw CG_FATAL( "Cards:handler" )
          << "Failed to parse the steering card at \"" << filename << "\"!\n"
          << err.what();
      }
    }

    void
    Handler::write( const Parameters* params, const std::string& filename )
    {
      try {
        auto writer = CardsHandlerFactory::get().build( extension( filename ) );
        writer->pack( params );
        return writer->write( filename );
      } catch ( const std::invalid_argument& err ) {
        throw CG_FATAL( "Cards:handler" )
          << "Failed to write the configuration to \"" << filename << "\"!\n"
          << err.what();
      }
    }
  }
}
