#include "CepGen/Core/Logger.h"

namespace CepGen
{
  Logger::Logger( std::ostream* os ) :
    level( Level::Information ), output( os )
  {}

  Logger::~Logger()
  {}

  void
  Logger::addExceptionRule( const std::string& rule )
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

  std::ostream&
  operator<<( std::ostream& os, const Logger::Level& lvl )
  {
    switch ( lvl ) {
      case Logger::Level::Nothing:
        return os << "None";
      case Logger::Level::Error:
        return os << "Errors";
      case Logger::Level::Warning:
        return os << "Warnings";
      case Logger::Level::Information:
        return os << "Infos";
      case Logger::Level::Debug:
        return os << "Debug";
      case Logger::Level::DebugInsideLoop:
        return os << "Debug (in loops)";
    }
    return os;
  }
}

