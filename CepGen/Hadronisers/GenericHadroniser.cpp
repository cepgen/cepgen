#include "CepGen/Hadronisers/GenericHadroniser.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  namespace Hadroniser
  {
    GenericHadroniser::GenericHadroniser( const char* name ) :
      name_( name )
    {}

    void
    GenericHadroniser::readStrings( const std::vector<std::string>& params ) {
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
  operator<<( std::ostream& os, const Hadroniser::GenericHadroniser& hadr )
  {
    return os << hadr.name().c_str();
  }

  std::ostream&
  operator<<( std::ostream& os, const Hadroniser::GenericHadroniser* hadr )
  {
    return os << hadr->name().c_str();
  }
}

