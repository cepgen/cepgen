#include "CepGen/Core/Logger.h"

namespace CepGen
{
  Logger::Logger( std::ostream* os ) :
    level( Level::information ), output( os )
  {}

  Logger::~Logger()
  {}

  void
  Logger::addExceptionRule( const std::string& rule )
  {
    allowed_exc_.emplace_back( rule );
  }

  bool
  Logger::passExceptionRule( const std::string& tmpl, const Level& lev ) const
  {
    if ( level >= lev )
      return true;
    if ( allowed_exc_.size() == 0 )
      return false;
    for ( const auto& rule : allowed_exc_ )
      if ( std::regex_match( tmpl, rule ) )
        return true;
    return false;
  }

  std::ostream&
  operator<<( std::ostream& os, const Logger::Level& lvl )
  {
    switch ( lvl ) {
      case Logger::Level::nothing:
        return os << "None";
      case Logger::Level::error:
        return os << "Errors";
      case Logger::Level::warning:
        return os << "Warnings";
      case Logger::Level::information:
        return os << "Infos";
      case Logger::Level::debug:
        return os << "Debug";
      case Logger::Level::debugInsideLoop:
        return os << "Debug (in loops)";
    }
    return os;
  }
}

