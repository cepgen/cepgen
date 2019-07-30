#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

namespace cepgen
{
  EventModifier::EventModifier( const ParametersList& plist, const std::string& name ) :
    name_( name ),
    seed_      ( plist.get<int>( "seed", -1ll ) ),
    max_trials_( plist.get<int>( "maxTrials", 1 ) )
  {
    CG_DEBUG( "EventModifier:init" )
      << "\"" << name_ << "\"-type event modifier built with:\n\t"
      << "* seed = " << seed_ << "\n\t"
      << "* maximum trials: " << max_trials_;
  }
}
