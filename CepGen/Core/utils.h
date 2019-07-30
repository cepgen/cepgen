#ifndef CepGen_Core_utils_h
#define CepGen_Core_utils_h

#include <string>
#include <vector>
#include <numeric>

namespace cepgen
{
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
  namespace utils
  {
    /// Add a trailing "s" when needed
    inline const char* s( size_t num ) { return ( num > 1 ) ? "s" : ""; }
    /// Add a trailing "s" when needed
    inline std::string s( const std::string& word, size_t num, bool show_number = false ) {
      return show_number
        ? Form( "%i %s%s", num, word.c_str(), ( num > 1 ) ? "s" : "" )
        : Form( "%s%s", word.c_str(), ( num > 1 ) ? "s" : "" );
    }
    /// Helper to print a vector
    template<class T> std::string repr( const std::vector<T>& vec, const std::string& sep = "," ) {
      return std::accumulate( std::next( vec.begin() ), vec.end(),
        std::to_string( *vec.begin() ), [&sep]( std::string str, T xv ) {
          return std::move( str )+sep+std::to_string( xv );
        } );
    }
    class ProgressBar
    {
      public:
        ProgressBar( size_t tot, size_t freq = 10 );
        void update( size_t iter ) const;

      private:
        static constexpr size_t BAR_LENGTH = 50;
        const std::string bar_pattern_;
        size_t total_, frequency_;
    };
  }
}

/// Provide a random number generated along a uniform distribution between 0 and 1
#define drand() (double)rand()/RAND_MAX

#endif
