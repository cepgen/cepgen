#include "CepGen/Core/utils.h"

#include <math.h>
#include <stdarg.h>  // For va_start, etc.

namespace CepGen
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
}
