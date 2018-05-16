#ifndef CepGen_Core_utils_h
#define CepGen_Core_utils_h

#include <string>
#include <stdarg.h>  // For va_start, etc.

namespace CepGen
{
  inline const char* s( unsigned short num ) { return ( num > 1 ) ? "s" : ""; }
}

/// Provide a random number generated along a uniform distribution between 0 and 1
//inline double drand() { srand (time(nullptr)); return static_cast<double>(rand())/RAND_MAX; }
#define drand() static_cast<double>( rand()/RAND_MAX )

/// Format a string using a printf style format descriptor.
inline std::string Form( const std::string fmt, ... ) {
  int size = ( (int)fmt.size() ) * 2 + 50;   // Use a rubric appropriate for your code
  std::string str;
  va_list ap;
  while ( true ) {     // Maximum two passes on a POSIX system...
    str.resize( size );
    va_start( ap, fmt );
    int n = vsnprintf( (char*)str.data(), size, fmt.c_str(), ap );
    va_end( ap );
    if ( n>-1 and n<size ) {  // Everything worked
      str.resize( n );
      return str;
    }
    if ( n>-1 )  // Needed size returned
      size = n + 1;   // For null char
    else size *= 2;      // Guess at a larger size (OS specific)
  }
  return str;
}


inline const char* yesno( const bool& test ) { return ( test ) ? "\033[32;1myes\033[0m" : "\033[31;1mno\033[0m"; }
//inline const char* boldify( const char* str ) { const std::string out = std::string( "\033[33;1m" ) + std::string( str ) + std::string( "\033[0m" ); return out.c_str(); }
inline std::string boldify( const std::string& str ) { return Form( "\033[1m%s\033[0m", str.c_str() ); }
inline std::string boldify( const char* str ) { return boldify( std::string( str ) ); }
inline std::string boldify( const double& dbl ) { return boldify( Form("%.2f", dbl ) ); }
inline std::string boldify( const int& i ) { return boldify( Form("% d", i ) ); }
inline std::string boldify( const unsigned int& ui ) { return boldify( Form("%d", ui ) ); }
inline std::string boldify( const unsigned long& ui ) { return boldify( Form("%lu", ui ) ); }
enum class Colour { gray = 30, red = 31, green = 32, yellow = 33, blue = 34, purple = 35 };
inline std::string colourise( const std::string& str, const Colour& col ) { return Form( "\033[%d%s\033[0m", (int)col, str.c_str() ); }

#endif
