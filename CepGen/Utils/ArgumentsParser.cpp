#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Core/Exception.h"

#include <sstream>
#include <algorithm>
#include <cstring>

namespace cepgen
{
  ArgumentsParser::ArgumentsParser( int argc, char* argv[] ) :
    help_str_( { { "help,h" } } ), config_str_( { { "cmd,c" } } ),
    help_req_( false )
  {
    command_name_ = argv[0];
    //--- first remove the program name
    std::vector<std::string> args_tmp;
    if ( argc > 1 ) {
      args_tmp.resize( argc-1 );
      std::copy( argv+1, argv+argc, args_tmp.begin() );
    }
    //--- then loop on user arguments to identify word -> value pairs
    for ( auto it_arg = args_tmp.begin(); it_arg != args_tmp.end(); ++it_arg ) {
      auto arg_val = utils::split( *it_arg, '=' ); // particular case for --arg=value
      //--- check if help message is requested
      for ( const auto& str : help_str_ )
        if ( arg_val.at( 0 ) == "--"+str.name.at( 0 )
          || ( str.name.size() > 1 && arg_val.at( 0 ) == "-"+str.name.at( 1 ) ) )
          help_req_ = true;
      //--- check if configuration word is requested
      for ( const auto& str : config_str_ )
        if ( arg_val.at( 0 ) == "--"+str.name.at( 0 )
          || ( str.name.size() > 1 && arg_val.at( 0 ) == "-"+str.name.at( 1 ) ) )
          extra_config_ = std::vector<std::string>( it_arg+1, args_tmp.end() );
      //--- parse arguments if word found after
      if ( arg_val.size() == 1 && arg_val.at( 0 )[0] == '-' && it_arg != std::prev( args_tmp.end() ) ) {
        const auto& word = *std::next( it_arg );
        if ( word[0] != '-' ) {
          arg_val.emplace_back( *std::next( it_arg ) );
          ++it_arg;
        }
      }
      args_.emplace_back( std::make_pair( arg_val.at( 0 ), arg_val.size() > 1 ? arg_val.at( 1 ) : "" ) );
    }
  }

  void
  ArgumentsParser::print_help() const
  {
    CG_INFO( "ArgumentsParser" ) << help_message();
  }

  void
  ArgumentsParser::dump() const
  {
    std::ostringstream os;
    os << "List of parameters retrieved from command-line:";
    for ( const auto& par : params_ )
      os
        << "\n\t[--" << par.name.at( 0 )
        << ( par.name.size() > 1 ? "|-"+par.name.at( 1 ) : "" )
        << ( par.optional ? ", optional" : "" )
        << "] = " << par.value;
    CG_INFO( "ArgumentsParser" ) << os.str();
  }

  ArgumentsParser&
  ArgumentsParser::parse()
  {
    if ( help_req_ ) {
      print_help();
      exit(0);
    }
    //--- loop over all parameters
    size_t i = 0;
    for ( auto& par : params_ ) {
      if ( par.name.empty() ) {
        //--- no argument name ; fetching by index
        if ( i >= args_.size() )
          throw CG_FATAL( "ArgumentsParser" )
            << help_message()
            << " Failed to retrieve required <arg" << i << ">.";
        par.value = !par.bool_variable ? args_.at( i ).second : "1";
      }
      else {
        // for each parameter, loop over arguments to find correspondance
        bool found_value = false;
        for ( const auto& arg : args_ ) {
          if ( arg.first == "--"+par.name.at( 0 )
            || ( par.name.size() > 1 && arg.first == "-"+par.name.at( 1 ) ) ) {
            par.value = arg.second;
            if ( par.bool_variable ) { // all particular cases for boolean arguments
              const auto word = utils::tolower( arg.second );
              if ( word.empty()
                || word == "on" || word != "off"
                || word == "yes" || word != "no"
                || word == "true" || word != "false" )
                par.value = "1"; // if the flag is set, enabled by default
            }
            found_value = true;
            ++i;
            break;
          }
        }
        if ( !found_value && args_.size() > i && args_.at( i ).first[0] != '-' )
          par.value = args_.at( i ).first;
        else if ( !found_value && !par.optional ) // no match
          throw CG_FATAL( "ArgumentsParser" )
            << help_message()
            << " The following parameter was not set: '" << par.name.at( 0 ) << "'.";
      }
      par.parse();
      CG_DEBUG( "ArgumentsParser" )
        << "Parameter '" << i << "|--" << par.name.at( 0 )
        << ( par.name.size() > 1 ? "|-"+par.name.at( 1 ) : "" ) << "'"
        << " has value '" << par.value << "'.";
      ++i;
    }
    return *this;
  }

  std::string
  ArgumentsParser::operator[]( std::string name ) const
  {
    for ( const auto& par : params_ ) {
      if ( "--"+par.name.at( 0 ) == name )
        return par.value;
      if ( par.name.size() > 1 && "-"+par.name.at( 1 ) == name )
        return par.value;
    }
    throw CG_FATAL( "ArgumentsParser" )
      << "The parameter \"" << name << "\" was not declared "
      << "in the arguments parser constructor!";
  }

  std::string
  ArgumentsParser::help_message() const
  {
    std::ostringstream oss;
    std::vector<std::pair<Parameter, size_t> > req_params, opt_params;
    oss << "Usage: " << command_name_;
    size_t i = 0;
    for ( const auto& par : params_ ) {
      if ( par.optional ) {
        opt_params.emplace_back( std::make_pair( par, i ) );
        oss << " [";
      }
      else {
        req_params.emplace_back( std::make_pair( par, i ) );
        oss << " ";
      }
      oss << ( !par.name.at( 0 ).empty() ? "--" : " <arg"+std::to_string( i )+">" )
        << par.name.at( 0 );
      if ( par.name.size() > 1 )
        oss << ( !par.name.at( 0 ).empty() ? "|" : "" ) << "-" << par.name.at( 1 );
      if ( par.optional )
        oss << "]";
      ++i;
    }
    if ( req_params.size() > 0 ) {
      oss << "\n    " << utils::s( "required argument", req_params.size(), false ) << ":";
      for ( const auto& par : req_params )
        oss << utils::format( "\n\t%s%-18s\t%-28s",
          ( par.first.name.size() > 1 ? "-"+par.first.name.at( 1 )+"|" : "" ).c_str(),
          ( !par.first.name.at( 0 ).empty() ? "--"+par.first.name.at( 0 ) : "<arg"+std::to_string( par.second )+">" ).c_str(),
          par.first.description.c_str() );
    }
    if ( opt_params.size() > 0 ) {
      oss << "\n    " << utils::s( "optional argument", opt_params.size(), false ) << ":";
      for ( const auto& par : opt_params )
        oss << utils::format( "\n\t%s%-18s\t%-28s\tdefault = '%s'",
          ( par.first.name.size() > 1 ? "-"+par.first.name.at( 1 )+"|" : "" ).c_str(),
          ( !par.first.name.at( 0 ).empty() ? "--"+par.first.name.at( 0 ) : "<arg"+std::to_string( par.second )+">" ).c_str(),
          par.first.description.c_str(), par.first.value.c_str() );
    }
    oss << std::endl;
    return oss.str();
  }

  //----- simple parameters

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::string* var, std::string default_value ) :
    name( utils::split( name, ',' ) ), description( description ),
    value( default_value ), optional( true ),
    str_variable( var ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, double* var, double default_value ) :
    name( utils::split( name, ',' ) ), description( description ),
    value( utils::format( "%g", default_value ) ), optional( true ),
    str_variable( nullptr ), float_variable( var ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, int* var, int default_value ) :
    name( utils::split( name, ',' ) ), description( description ),
    value( utils::format( "%+i", default_value ) ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( var ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, unsigned int* var, unsigned int default_value ) :
    name( utils::split( name, ',' ) ), description( description ),
    value( std::to_string( default_value ) ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( var ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, bool* var, bool default_value ) :
    name( utils::split( name, ',' ) ), description( description ),
    value( utils::format( "%d", default_value ) ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( var ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  //----- vector of parameters

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<std::string>* var, std::vector<std::string> default_value ) :
    name( utils::split( name, ',' ) ), description( description ),
    value( utils::merge( default_value, "," ) ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( var ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<int>* var, std::vector<int> default_value ) :
    name( utils::split( name, ',' ) ), description( description ),
    value( "" ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( var ), vec_float_variable( nullptr )
  {
    std::string sep;
    for ( const auto& val : default_value )
      value += sep+utils::format( "%d", val ), sep = ",";
  }

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<double>* var, std::vector<double> default_value ) :
    name( utils::split( name, ',' ) ), description( description ),
    value( "" ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( var )
  {
    std::string sep;
    for ( const auto& val : default_value )
      value += sep+utils::format( "%e", val ), sep = ",";
  }

  ArgumentsParser::Parameter&
  ArgumentsParser::Parameter::parse()
  {
    CG_DEBUG( "ArgumentsParser:Parameter:parse" )
      << "Parsing argument " << name << ".";
    if ( str_variable != nullptr ) {
      *str_variable = value;
      return *this;
    }
    if ( float_variable != nullptr )
      try {
        *float_variable = std::stod( value );
        return *this;
      } catch ( const std::invalid_argument& ) {
        throw CG_FATAL( "ArgumentsParser:Parameter:parse" )
          << "Failed to parse variable '" << name << "' as float!";
      }
    if ( int_variable != nullptr )
      try {
        *int_variable = std::stoi( value );
        return *this;
      } catch ( const std::invalid_argument& ) {
        throw CG_FATAL( "ArgumentsParser:Parameter:parse" )
          << "Failed to parse variable '" << name << "' as integer!";
      }
    if ( uint_variable != nullptr )
      try {
        *uint_variable = std::stoi( value );
        return *this;
      } catch ( const std::invalid_argument& ) {
        throw CG_FATAL( "ArgumentsParser:Parameter:parse" )
          << "Failed to parse variable '" << name << "' as unsigned integer!";
      }
    if ( bool_variable != nullptr ) {
      try {
        *bool_variable = ( std::stoi( value ) != 0 );
        return *this;
      } catch ( const std::invalid_argument& ) {
        *bool_variable = ( strcasecmp( "true", value.c_str() ) == 0
                        || strcasecmp( "yes", value.c_str() ) == 0
                        || strcasecmp( "on", value.c_str() ) == 0 )
                        && strcasecmp( "false", value.c_str() ) != 0
                        && strcasecmp( "no", value.c_str() ) != 0
                        && strcasecmp( "off", value.c_str() ) != 0;
      }
    }
    if ( vec_str_variable != nullptr ) {
      *vec_str_variable = utils::split( value, ',' );
      return *this;
    }
    if ( vec_int_variable != nullptr ) {
      vec_int_variable->clear();
      const auto buf = utils::split( value, ',' );
      std::transform( buf.begin(), buf.end(), std::back_inserter( *vec_int_variable ), []( const std::string& str ) { return std::stoi( str ); } );
      return *this;
    }
    if ( vec_float_variable != nullptr ) {
      vec_float_variable->clear();
      const auto buf = utils::split( value, ',' );
      std::transform( buf.begin(), buf.end(), std::back_inserter( *vec_float_variable ), []( const std::string& str ) { return std::stod( str ); } );
      return *this;
    }
    throw CG_FATAL( "Parameter" )
      << "Failed to parse parameter \"" << name.at( 0 ) << "\"!";
  }

  std::ostream&
  operator<<( std::ostream& os, const ArgumentsParser::Parameter& par )
  {
    return os << "Parameter{"
      << "--" << par.name.at( 0 )
      << ( par.name.size() > 1 ? ",-"+par.name.at( 1 ) : "" )
      << ( !par.description.empty() ? ","+par.description : "" )
      << ",val=" << par.value
      << ",opt:" << std::boolalpha << par.optional
      << "}";
  }
}

