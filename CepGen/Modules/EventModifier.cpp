#include "CepGen/Modules/EventModifier.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

namespace cepgen
{
  EventModifier::EventModifier( const ParametersList& plist ) :
    name_( plist.get<std::string>( ParametersList::MODULE_NAME, "<invalid>" ) ),
    seed_( plist.get<int>( "seed", -1ll ) ),
    max_trials_( plist.get<int>( "maxTrials", 1 ) )
  {
    CG_DEBUG( "EventModifier:init" )
      << "\"" << name_ << "\"-type event modifier built with:\n\t"
      << "* seed = " << seed_ << "\n\t"
      << "* maximum trials: " << max_trials_;
  }

  void
  EventModifier::readStrings( const std::vector<std::string>& params )
  {
    if ( params.empty() )
      return;
    std::ostringstream os;
    for ( const auto& p : params ) {
      readString( p );
      os << "\n\t  '" << p << "'";
    }
    CG_DEBUG( "EventModifier:configure" )
      << "Feeding \"" << name_ << "\" event modifier algorithm with:"
      << os.str();
  }
}
