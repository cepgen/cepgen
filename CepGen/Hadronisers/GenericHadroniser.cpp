#include "CepGen/Hadronisers/GenericHadroniser.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

namespace cepgen
{
  namespace hadroniser
  {
    GenericHadroniser::GenericHadroniser( const char* name, const ParametersList& plist ) :
      name_( name ),
      seed_      ( plist.get<int>( "seed", -1ll ) ),
      max_trials_( plist.get<int>( "maxTrials", 1 ) )
    {
      CG_DEBUG( "Hadroniser:init" )
        << "\"" << name_ << "\"-type hadroniser built with:\n\t"
        << "* seed = " << seed_ << "\n\t"
        << "* maximum trials: " << max_trials_;
    }

    void
    GenericHadroniser::readStrings( const std::vector<std::string>& params ) {
      if ( params.empty() )
        return;
      std::ostringstream os;
      for ( const auto& p : params ) {
        readString( p );
        os << "\n\t  '" << p << "'";
      }
      CG_DEBUG( "Hadroniser:configure" )
        << "Feeding \"" << name_ << "\" hadroniser with:"
        << os.str();
    }

    std::string
    GenericHadroniser::name() const
    {
      return name_;
    }
  }

  std::ostream&
  operator<<( std::ostream& os, const hadroniser::GenericHadroniser& hadr )
  {
    return os << hadr.name().c_str();
  }

  std::ostream&
  operator<<( std::ostream& os, const hadroniser::GenericHadroniser* hadr )
  {
    return os << hadr->name().c_str();
  }
}
