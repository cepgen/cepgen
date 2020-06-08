#include "CepGen/Utils/String.h"

#include <iomanip>
#include <iterator>
#include <algorithm>
#include <sstream>

#include <math.h>
#include <stdarg.h>  // For va_start, etc.

namespace cepgen
{
  namespace utils
  {
    std::string
    format( const std::string fmt, ... )
    {
      int size = ( (int)fmt.size() ) * 2 + 50;
      std::string str;
      va_list ap;
      while ( true ) {
        //--- maximum two passes on a POSIX system...
        str.resize( size );
        va_start( ap, fmt );
        int n = vsnprintf( (char*)str.data(), size, fmt.c_str(), ap );
        va_end( ap );
        //--- check if everything worked
        if ( n > -1 && n < size ) {
          str.resize( n );
          return str;
        }
        size = ( n > -1 ) ? n+1 : size*2;
      }
      return str;
    }

    std::string
    yesno( bool test )
    {
      return test
        ? "\033[32;1myes\033[0m"
        : "\033[31;1mno\033[0m";
    }

    /// String implementation of the boldification procedure
    template<> std::string
    boldify<std::string>( std::string str )
    {
      return format( "\033[1m%s\033[0m", str.c_str() );
    }

    /// C-style character array implementation of the boldification procedure
    template<> std::string
    boldify<const char*>( const char* str )
    {
      return boldify( std::string( str ) );
    }

    /// Unsigned long integer implementation of the boldification procedure
    template<> std::string
    boldify( unsigned long ui )
    {
      return boldify( std::to_string( ui ) );
    }

    std::string
    colourise( const std::string& str, const Colour& col )
    {
      return format( "\033[0;%dm%s\033[0m", (int)col, str.c_str() );
    }

    size_t
    replace_all( std::string& str, const std::string& from, const std::string& to )
    {
      size_t count = 0, pos = 0;
      while ( ( pos = str.find( from, pos ) ) != std::string::npos ) {
        str.replace( pos, from.length(), to );
        pos += to.length();
        ++count;
      }
      return count;
    }

    std::vector<std::string>
    split( const std::string& str, char delim )
    {
      std::vector<std::string> out;
      std::string token;
      std::istringstream iss( str );
      while ( std::getline( iss, token, delim ) )
        out.emplace_back( token );
      return out;
    }

    std::string
    merge( const std::vector<std::string>& vec, const std::string& delim )
    {
      if ( vec.empty() )
        return std::string();
      if ( vec.size() == 1 )
        return vec.at( 0 );
      std::ostringstream oss;
      std::copy( vec.begin(), std::prev( vec.end() ), std::ostream_iterator<std::string>( oss, delim.c_str() ) );
      return oss.str()+*vec.rbegin();
    }

    std::string
    toupper( const std::string& str )
    {
      std::string out; out.resize( str.size() );
      std::transform( str.begin(), str.end(), out.begin(), ::toupper );
      return out;
    }

    std::string
    tolower( const std::string& str )
    {
      std::string out; out.resize( str.size() );
      std::transform( str.begin(), str.end(), out.begin(), ::tolower );
      return out;
    }

    std::string
    environ( const std::string& env, const std::string& def )
    {
      const auto out = std::getenv( env.c_str() );
      if ( !out )
        return def;
      return std::string( out );
    }
  }
}
