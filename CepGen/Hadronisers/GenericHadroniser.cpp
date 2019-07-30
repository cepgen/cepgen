#include "CepGen/Hadronisers/GenericHadroniser.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

namespace cepgen
{
  namespace hadr
  {
    GenericHadroniser::GenericHadroniser( const ParametersList& plist, const std::string& name ) :
      EventModifier( plist, name ),
      remn_fragm_( plist.get<bool>( "remnantsFragmentation", true ) )
    {}

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
  }

  std::ostream&
  operator<<( std::ostream& os, const hadr::GenericHadroniser& hadr )
  {
    return os << hadr.name();
  }

  std::ostream&
  operator<<( std::ostream& os, const hadr::GenericHadroniser* hadr )
  {
    return os << hadr->name();
  }
}
