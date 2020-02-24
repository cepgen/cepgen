#include "CepGen/Cards/Handler.h"
#include "CepGen/Modules/CardsHandlerFactory.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  namespace card
  {
    /// Command line parser
    class CommandLineHandler : public Handler
    {
      public:
        explicit CommandLineHandler( const ParametersList& );

      private:
        const std::vector<std::string> argv_;
    };

    CommandLineHandler::CommandLineHandler( const ParametersList& params ) :
      argv_( params.get<std::vector<std::string> >( "args" ) )
    {
      for ( const auto& arg : argv_ )
        CG_INFO("") << ">>> " << arg;
    }
  }
}

REGISTER_CARD_HANDLER( "cmd", CommandLineHandler )
