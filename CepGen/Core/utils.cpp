#include "CepGen/Core/utils.h"

#include <iomanip>
#include <iterator>
#include <algorithm>

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
    std::ostringstream oss;
    std::copy( vec.begin(), std::prev( vec.end() ), std::ostream_iterator<std::string>( oss, delim.c_str() ) );
    return oss.str()+*vec.rbegin();
  }

  namespace utils
  {
    ProgressBar::ProgressBar( size_t tot, size_t freq ) :
      bar_pattern_( BAR_LENGTH, '=' ), total_( tot ), frequency_( freq )
    {}

    void ProgressBar::update( size_t iter ) const
    {
      const size_t percent = iter*100./total_;
      if ( percent % frequency_ == 0 || iter == total_ ) {
        int lpad = int( percent/100. * BAR_LENGTH );
        int rpad = BAR_LENGTH-lpad;
        fprintf( stderr, "\r%3zu%% [%.*s%*s]", percent, lpad, bar_pattern_.c_str(), rpad, "" );
        fflush( stderr );
        if ( iter == total_ )
          fprintf( stderr, "\n" );
      }
    }
  }
}
