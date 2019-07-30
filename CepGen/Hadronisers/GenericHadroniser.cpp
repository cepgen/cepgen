#include "CepGen/Hadronisers/GenericHadroniser.h"
#include "CepGen/Core/ParametersList.h"

namespace cepgen
{
  namespace hadr
  {
    GenericHadroniser::GenericHadroniser( const ParametersList& plist, const std::string& name ) :
      EventModifier( plist, name ),
      remn_fragm_( plist.get<bool>( "remnantsFragmentation", true ) )
    {}
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
