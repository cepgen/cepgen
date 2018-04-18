#include "CepGen/Core/Logger.h"

namespace CepGen
{
  void
  Logger::addExceptionRule( const char* rule )
  {
    allowed_exc_.emplace_back( rule );
  }

  bool
  Logger::passExceptionRule( const std::string& tmpl ) const
  {
    for ( const auto& rule : allowed_exc_ )
      if ( std::regex_match( tmpl, rule ) )
        return true;
    return false;
  }
}

