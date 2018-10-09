#ifndef CepGen_Core_utils_h
#define CepGen_Core_utils_h

#include <string>

namespace cepgen
{
  /// Add a closing "s" when needed
  inline const char* s( unsigned short num ) { return ( num > 1 ) ? "s" : ""; }
  /// Format a string using a printf style format descriptor.
  std::string Form( const std::string fmt, ... );
  /// Human-readable boolean printout
  inline const char* yesno( const bool& test ) { return ( test ) ? "\033[32;1myes\033[0m" : "\033[31;1mno\033[0m"; }
  //inline const char* boldify( const char* str ) { const std::string out = std::string( "\033[33;1m" ) + std::string( str ) + std::string( "\033[0m" ); return out.c_str(); }
  /// Boldify a string for TTY-type output streams
  inline std::string boldify( const std::string& str ) { return Form( "\033[1m%s\033[0m", str.c_str() ); }
  /// Boldify a string for TTY-type output streams
  inline std::string boldify( const char* str ) { return boldify( std::string( str ) ); }
  /// Boldify a double floating point number for TTY-type output streams
  inline std::string boldify( const double& dbl ) { return boldify( Form("%.2f", dbl ) ); }
  /// Boldify an integer for TTY-type output streams
  inline std::string boldify( const int& i ) { return boldify( Form("% d", i ) ); }
  /// Boldify an unsigned integer for TTY-type output streams
  inline std::string boldify( const unsigned int& ui ) { return boldify( Form("%d", ui ) ); }
  /// Boldify an unsigned long integer for TTY-type output streams
  inline std::string boldify( const unsigned long& ui ) { return boldify( Form("%lu", ui ) ); }
  /// TTY-type enumeration of colours
  enum class Colour { gray = 30, red = 31, green = 32, yellow = 33, blue = 34, purple = 35 };
  /// Colourise a string for TTY-type output streams
  inline std::string colourise( const std::string& str, const Colour& col ) { return Form( "\033[%d%s\033[0m", (int)col, str.c_str() ); }
  /// Replace all occurences of a text by another
  size_t replace_all( std::string& str, const std::string& from, const std::string& to );
}

/// Provide a random number generated along a uniform distribution between 0 and 1
#define drand() static_cast<double>( rand()/RAND_MAX )

#endif
