#include "CepGen/Utils/String.h"

#include <iomanip>
#include <iterator>
#include <algorithm>
#include <sstream>

#include <math.h>
#include <stdarg.h>  // For va_start, etc.

namespace cepgen
{
  std::string
  Form( const std::string fmt, ... )
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
}
